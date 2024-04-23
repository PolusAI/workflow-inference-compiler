# Compute API

As previously mentioned, one of the beautiful things about the declarative approach to workflows is that we can execute workflows on massive machines just as easily as executing workflows on a local laptop. Concretely, merely changing `--run_local` to `--run_compute`, we can execute the exact same workflow on the NCATS HPC cluster! That's it! Absolutely no modifications to the workflow itself are necessary!

Note that if your workflow inputs file references files on your local machine, then
1. those files will need to be manually transferred to the new machine (this could perhaps be automated, but maybe you have ~1TB of data...)
2. if you are using absolute file paths, the paths will need to be updated w.r.t. the new machine. (So use relative paths ;) )

## Authentication Access Token

When using `--run_compute` you will also need to use `--compute_access_token $ACCESS_TOKEN`. Unfortunately, there is currently no programmatic way of obtaining the access token via an API from the command line. You will need to manually perform the following steps:

* Go to https://compute.scb-ncats.io/explorer/
* Click the green Authorize button. You will be taken to the NIH login page.
* Enter your NIH username and password, then
* Authenticate (using Microsoft Authenticator)
* You will be returned to https://compute.scb-ncats.io/explorer/
* Click Close (NOT Logout!)
* Scroll down to [HealthController](https://compute.scb-ncats.io/explorer/#/HealthController/HealthController.ping)
* Click Try It Out and then click Execute.
* You should see a massive hash string after "authorization: bearer".
* Copy the massive hash string. Be careful not to copy any other characters (like blank spaces).
* In a bash terminal, create the environment variable `export ACCESS_TOKEN=...` (where ... means paste in the hash string)

Unfortunately, the access token currently expires after about an hour, so you will need to repeat these steps periodically.

## Workflow Status

After submitting a workflow, users can check on its status either using https://compute.scb-ncats.io/explorer/ or by directly logging into dali.ncats.nih.gov and using the `squeue` command. Currently, the workflows are executed under the placeholder svc-polus user.

```
ssh <username>@dali.ncats.nih.gov
```
NOTE: The above server is behind the NIH VPN; you must enable the VPN to access it.

```
watch -n5 squeue -u svc-polus
```

After a ~2-3 minute initial delay you should see nodes starting, running, and finishing, corresponding to the individual steps in the workflow. The output files and logs are currently stored under /project/labshare-compute/