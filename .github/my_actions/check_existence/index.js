// See https://docs.github.com/en/actions/creating-actions/creating-a-javascript-action
// NOTE: Every time you modify this file, you need to run
// `ncc build index.js && git add -f dist/index.js index.js package.json package-lock.json`
// Install ncc using `npm install -g @vercel/ncc`
// You do NOT need to git add node_modules/*

const core = require('@actions/core');
const github = require('@actions/github');
import fetch from "node-fetch";

try {
  const repository = core.getInput('repository');
  const sender_repo_owner = core.getInput('sender_repo_owner');
  const sender_repo_ref = core.getInput('sender_repo_ref');
  const default_owner = core.getInput('default_owner');
  const default_branch = core.getInput('default_branch');
  const access_token = core.getInput('access_token');

  if (!access_token) {
    console.log("Error! access_token is not defined! (or expired)");
  }

  const url_branches = "https://api.github.com/repos/" + sender_repo_owner + "/" + repository + "/branches";
  console.log(`url_branches: ${url_branches}`);
  console.log(`access_token: ${access_token}`);


  // https://docs.github.com/en/rest/overview/resources-in-the-rest-api?apiVersion=2022-11-28#rate-limiting
  // NOTE: Even though this does not require authentication, unauthenticated requests
  // are limited to 60 per hour, whereas authenticated requests are 5000 per hour.
  const response = await fetch(url_branches, {
    method: "GET",
    headers: {
      //      "Accept": "application/vnd.github+json",
      //      "X-GitHub-Api-Version": "2022-11-28",
            'Authorization': `Bearer ${access_token}`
          }
  });
  const branches = await response.json();

  var exists = false;
  if (branches.message == "Not Found") {
    console.log(`${url_branches} not found.`);
    const url_default = "https://github.com/" + default_owner + "/" + repository;
    console.log(`This probably means you did not fork the repo: ${url_default}`);
    console.log(`Falling back to ${url_default}`);
  } else {
    const branches_str = JSON.stringify(branches, undefined, 2)
    console.log(`response.json(): ${branches_str}`);

    const equals_name = (branch) => branch.name === sender_repo_ref;
    exists = branches.some(equals_name);
  }
  const owner = exists ? sender_repo_owner : default_owner;
  const ref = exists ? sender_repo_ref : default_branch;

  core.setOutput("repository", repository);
  core.setOutput("owner", owner);
  core.setOutput("ref", ref);

  // Get the JSON webhook payload for the event that triggered the workflow
  const payload = JSON.stringify(github.context.payload, undefined, 2)
  console.log(`The event payload: ${payload}`);
} catch (error) {
  core.setFailed(error.message);
}
