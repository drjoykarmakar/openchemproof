from types import SimpleNamespace
from unittest.mock import AsyncMock

import pytest

from benchmarks.chem_artifact_bench.eval_runner import evaluate_one
from benchmarks.chem_artifact_bench.generate_bench import generate_cases


def completion(content=None, tool_calls=None):
    openai = pytest.importorskip("openai")
    return openai.types.chat.ChatCompletion.model_validate(
        {
            "id": "mock-response",
            "model": "mock-model",
            "created": 0,
            "object": "chat.completion",
            "choices": [
                {
                    "index": 0,
                    "finish_reason": "stop",
                    "message": {"role": "assistant", "content": content, "tool_calls": tool_calls},
                }
            ],
        }
    )


async def test_mocked_baseline_has_no_tools_or_gold_labels():
    create = AsyncMock(return_value=completion('{"status":"INCONCLUSIVE","rationale":"missing evidence"}'))
    api = SimpleNamespace(chat=SimpleNamespace(completions=SimpleNamespace(create=create)))
    case = generate_cases()[0]
    row = await evaluate_one(api, "mock-model", case, "baseline")
    assert row["status"] == "INCONCLUSIVE" and "error" not in row
    arguments = create.call_args.kwargs
    assert "tools" not in arguments
    prompt = arguments["messages"][1]["content"]
    assert case.category not in prompt and "expected_status" not in prompt


async def test_mocked_live_loop_really_calls_mcp():
    case = generate_cases()[0]
    import json

    call = {
        "id": "call-1",
        "type": "function",
        "function": {
            "name": "run_full_verification_pipeline",
            "arguments": json.dumps({"record": case.record.model_dump(mode="json")}),
        },
    }
    create = AsyncMock(
        side_effect=[
            completion(tool_calls=[call]),
            completion('{"status":"ARTIFACT_DETECTED","rationale":"tool evidence"}'),
        ]
    )
    api = SimpleNamespace(chat=SimpleNamespace(completions=SimpleNamespace(create=create)))
    row = await evaluate_one(api, "mock-model", case, "tool_augmented")
    assert "error" not in row
    tools = [item for item in row["trace"] if "tool" in item]
    assert len(tools) == 1
    assert "audit_sha256" in json.loads(tools[0]["result"])


async def test_malformed_model_answer_is_recorded():
    create = AsyncMock(return_value=completion("not json"))
    api = SimpleNamespace(chat=SimpleNamespace(completions=SimpleNamespace(create=create)))
    row = await evaluate_one(api, "mock-model", generate_cases()[0], "baseline")
    assert row["status"] == "INCONCLUSIVE" and row["error"]
