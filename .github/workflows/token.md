1. Automatic authentication using "installation access token" of Github App needs
inputs of the App ID and the App private key. Here's an example of automating authentication
with Github App:
https://docs.github.com/en/issues/planning-and-tracking-with-projects/automating-your-project/automating-projects-using-actions
2. The App ID and App private key are stored as action secrets of the repository.
Collaborators of public fork of a personal account won't have access to these secrets, even
though they have write permissions of the fork. The secrets are accessible for users with
"maint"
3. This third-party Github Action repo is used in the above example of the official Github
doc. To prevent future updates of the Action repo to change behaviors and posing security
risks, explicity use the commit hash of the most recent stable release (v1.8.0) to freeze
the used code.
* v1.8.0 of the third-party Action repo:
  https://github.com/tibdex/github-app-token/tree/releases/v1.8.0