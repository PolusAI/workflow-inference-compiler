// See https://docs.github.com/en/actions/creating-actions/creating-a-javascript-action
// NOTE: Every time you modify this file, you need to run
// `ncc build index.js && git add -f dist/index.js index.js package.json package-lock.json`
// You do NOT need to git add node_modules/*

const core = require('@actions/core');
const github = require('@actions/github');
import fetch from "node-fetch";

try {
  const github_repository = core.getInput('github_repository');
  const event_type = core.getInput('event_type');
  const sender_repo = core.getInput('sender_repo');
  const operating_system = core.getInput('operating_system');
  const commit_message = core.getInput('commit_message');
  const mm_workflows_ref = core.getInput('mm_workflows_ref');
  const workflow_success = core.getInput('workflow_success');
  const access_token = core.getInput('access_token');

  if (!access_token) {
    console.log("Error! access_token is not defined! (or expired)");
  }

  const url_dispatches = "https://api.github.com/repos/" + sender_repo + "/dispatches";
  console.log(`url_branches: ${url_dispatches}`);
  console.log(`access_token: ${access_token}`);

  const run_name = commit_message + " - " + mm_workflows_ref + " - " + operating_system;

  if (event_type == "repository_dispatch") {
    console.log(`Replying to sender_repo: ${sender_repo}`);
    const response = await fetch(url_dispatches, {
      method: "POST",
      body: JSON.stringify({
        event_type: "ci_status",
        client_payload: {
          repository: github_repository,
          success: workflow_success,
          run_name: run_name,
        },
      }),
      headers: {
  //      "Accept": "application/vnd.github+json",
  //      "X-GitHub-Api-Version": "2022-11-28",
        'Authorization': `Bearer ${access_token}`
      }
    });
    const response_str = await response.text();
    console.log(`response: ${response_str}`);
    core.setOutput("response", response_str);
  }

  // Get the JSON webhook payload for the event that triggered the workflow
  const payload = JSON.stringify(github.context.payload, undefined, 2)
  console.log(`The event payload: ${payload}`);
} catch (error) {
  core.setFailed(error.message);
}
