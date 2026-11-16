import asyncio
import os
import socket
import subprocess
import sys

import httpx2
from fastmcp import Client
from fastmcp.client.transports import StdioTransport, StreamableHttpTransport


async def test_stdio_entrypoint(ethanol, tmp_path):
    transport = StdioTransport(
        command=sys.executable,
        args=["-m", "openchemproof.mcp_server.server"],
        env=dict(os.environ),
        log_file=tmp_path / "stdio.log",
    )
    async with Client(transport) as client:
        assert len(await client.list_tools()) == 5
        result = await client.call_tool("verify_structure_identity", {"record": ethanol.model_dump(mode="json")})
        assert result.structured_content["status"] == "INCONCLUSIVE"


async def test_streamable_http_entrypoint(ethanol, tmp_path):
    with socket.socket() as listener:
        listener.bind(("127.0.0.1", 0))
        port = listener.getsockname()[1]
    with (tmp_path / "http.log").open("w") as log:
        process = subprocess.Popen(
            [sys.executable, "-m", "openchemproof.mcp_server.server", "--transport", "http", "--port", str(port)],
            stdout=log,
            stderr=log,
        )
        try:
            for _ in range(100):
                assert process.poll() is None, "HTTP server exited unexpectedly"
                try:
                    with socket.create_connection(("127.0.0.1", port), timeout=0.1):
                        break
                except OSError:
                    await asyncio.sleep(0.1)
            # The test connects only to its own loopback server, independent of host proxy settings.
            transport = StreamableHttpTransport(
                f"http://127.0.0.1:{port}/mcp",
                httpx_client_factory=lambda **kwargs: httpx2.AsyncClient(**kwargs, trust_env=False),
            )
            async with Client(transport) as client:
                assert len(await client.list_resources()) == 2
                result = await client.call_tool(
                    "run_full_verification_pipeline", {"record": ethanol.model_dump(mode="json")}
                )
                assert result.structured_content["record_id"] == ethanol.record_id
        finally:
            process.terminate()
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5)
