A commit within a single repository is an atomic transaction of changes.
To define an atomic "transaction" of commits across repositories,
we can simply create branches in all of the repos with the exact same name.
Then, on push or pull request to any individual repo, we query all of the
repos for a branch with that exact same name, and run all of the individual CIs
using those particular branches. If a branch isn't found in a given repo,
we fallback to the default branch on the upstream repo (i.e. PolusAI & master/main)
(So you should always use unique branch names; do not use 'testing', etc across repos.)

We are linking the CIs using asynchronous API calls, so they will show up in
your actions histories (plural!) as separate runs. Moreover, actions initiated
by repository_dispatch or workflow_dispatch are NOT affiliated with any branch
(How could they be? Anyone or anything can call the API...)
so the github website will NOT display the dispatching branch name in blue.
In particular, for PRs,

<span style="color:red"> A GREEN CHECK MARK FOR THE INITIAL DISPATCH WORKFLOW DOES NOT MEAN ANY OF THE OTHER ACTIONS PASSED! </span>

<span style="color:red"> YOU MUST MANUALLY EXAMINE THE ACTIONS HISTORIES OF *ALL* LINKED REPOS! </span>

We have very carefully passed along the commit message, repository, and branch
names so users can easily see which actions correspond to a given
"transaction" / coordinated set of branches across repositories.
Note that since the API calls are asynchronous, there are unavoidable race conditions.
In particular, when merging PRs in different repos as part of a "transaction",
you have about 10 seconds to simultaneously merge all PRs into master.
Otherwise, you may simply need to manually re-run the actions.

Another very important point about PRs is that the authentication currently uses
Personal Access Tokens (see below), so the API calls will all fail on PRs.
However, if all of the branches are up-to-date with all of the upstream master
branches, then if the CIs all pass on all of the forks, then we don't necessarily
need to see a green check mark on the PRs. This approach will require
frequent rebasing, but that's a good thing.

The purpose of this file is to contain (almost) all of the logic necessary
for linking CIs across repos. Once the appropriate forks and branches are
determined here, you can just pass them as inputs into your existing github
actions workflows, and then add a few steps at the end to return the result.

First note that we can only call the default branch (master/main) through repository_dispatch.
Note that is is not sufficient to checkout master/main and then `git switch`
to another branch; the actions workflow files themselves may be different
between branches! So after we find the forks and branches here, we can use
workflow_dispatch to directly checkout the correct code and run the CIs on
the correct forks/branches.
https://github.com/orgs/community/discussions/24657#discussioncomment-3244904
https://docs.github.com/en/rest/actions/workflows?apiVersion=2022-11-28#create-a-workflow-dispatch-event

Each type of webhook events has its own set of body parameters.
'repository_dispatch' requires 'event_type' and uses 'client_payload' to pass user data
https://docs.github.com/en/webhooks-and-events/webhooks/webhook-events-and-payloads#repository_dispatch
'workflow_dispatch' requires 'ref' and uses 'inputs' to pass user data.
https://docs.github.com/en/webhooks-and-events/webhooks/webhook-events-and-payloads#workflow_dispatch

Finally, in general rather than using "$GITHUB_*", prefer the equivalent "${{ github.* }}".
https://docs.github.com/en/actions/learn-github-actions/variables#default-environment-variables
https://docs.github.com/en/actions/learn-github-actions/contexts#github-context
