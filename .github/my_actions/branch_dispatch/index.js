// See https://docs.github.com/en/actions/creating-actions/creating-a-javascript-action
// NOTE: Every time you modify this file, you need to run
// `ncc build index.js && git add -f dist/index.js index.js package.json package-lock.json`
// You do NOT need to git add node_modules/*

const core = require('@actions/core');
const github = require('@actions/github');
import fetch from "node-fetch";

try {
  const repository = core.getInput('repository');
  const workflow_yml = core.getInput('workflow_yml');
  const sender_repo = core.getInput('sender_repo');
  const sender_repo_owner = core.getInput('sender_repo_owner');
  const wic_owner = core.getInput('wic_owner');
  const wic_ref = core.getInput('wic_ref');
  const event_type = core.getInput('event_type');
  const commit_message = core.getInput('commit_message');
  const mm_workflows_owner = core.getInput('mm_workflows_owner');
  const mm_workflows_ref = core.getInput('mm_workflows_ref');
  const access_token = core.getInput('access_token');

  if (!access_token) {
    console.log("Error! secrets.ACCESS_TOKEN is not defined! (or expired)");
  }

  // Use base repository owner, otherwise permission errors: Resource not accessible by integration.
  const repo_owner = event_type == "pull_request_target" ? github.context.payload.repository.owner.login : sender_repo_owner;
  const url_dispatches = "https://api.github.com/repos/" + repo_owner + "/" + repository + "/actions/workflows/" + workflow_yml + "/dispatches";
  console.log(`url_branches: ${url_dispatches}`);
  console.log(`access_token: ${access_token}`);

  const response = await fetch(url_dispatches, {
    method: "POST",
    body: JSON.stringify({
      ref: wic_ref,
      inputs: {
        event_type: event_type,
        commit_message: commit_message,
        sender_repo: sender_repo,
        "wic_owner": wic_owner,
        "wic_ref": wic_ref,
        "mm-workflows_owner": mm_workflows_owner,
        "mm-workflows_ref": mm_workflows_ref,
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
  core.setOutput("respose", response_str);

  // Get the JSON webhook payload for the event that triggered the workflow
  const payload = JSON.stringify(github.context.payload, undefined, 2)
  console.log(`The event payload: ${payload}`);
} catch (error) {
  core.setFailed(error.message);
}
