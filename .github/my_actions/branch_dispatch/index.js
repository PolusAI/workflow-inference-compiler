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
  const dispatch_ref = core.getInput('dispatch_ref');
  const wic_owner = core.getInput('wic_owner');
  const wic_ref = core.getInput('wic_ref');
  const event_type = core.getInput('event_type');
  const commit_message = core.getInput('commit_message');
  const mm_workflows_owner = core.getInput('mm_workflows_owner');
  const mm_workflows_ref = core.getInput('mm_workflows_ref');
  const access_token = core.getInput('access_token');

  if (!access_token) {
    console.log("Error! access_token is not defined! (or expired)");
  }

  // Use base repository owner, otherwise permission errors: Resource not accessible by integration.
  const url_dispatches = "https://api.github.com/repos/" + github.context.repo.owner + "/" + repository + "/actions/workflows/" + workflow_yml + "/dispatches";
  console.log(`url_branches: ${url_dispatches}`);
  console.log(`access_token: ${access_token}`);

  // Note: In PRs, the branch to which the dispatch is sent to is not necessarily the same
  // as the branch to be checked out. Especially for the case of upstream being PolusAI, the
  // dispatch should be sent to the 'master' branch while the code to be checked out should be
  // the feature branch of the user's fork. When seeing error message: "message":"No ref found for: ..."
  // the problem is in fact that (with the second form of authentication) the username is wrong
  // (i.e. username in dispatch URL being the upstream).
  const response = await fetch(url_dispatches, {
    method: "POST",
    body: JSON.stringify({
      ref: dispatch_ref,
      inputs: {
        "event_type": event_type,
        "commit_message": commit_message,
        "sender_repo": sender_repo,
        "sender_repo_owner": sender_repo_owner,
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
  core.setOutput("response", response_str);

  // Get the JSON webhook payload for the event that triggered the workflow
  const payload = JSON.stringify(github.context.payload, undefined, 2)
  console.log(`The event payload: ${payload}`);
} catch (error) {
  core.setFailed(error.message);
}
