# Running on HPC (SLURM) and AWS

Scripts and settings below were validated locally (see "What was verified").

Needs outbound HTTPS to rest.ensembl.org and GitHub (MHCflurry weights); many HPC compute nodes have no internet, so use a login-node session or pre-stage the weights.

## HPC (SLURM + Apptainer)

Uses the `slurm` profile in `nextflow.config` (Singularity/Apptainer
containers). Edit `process.queue` for your cluster, then:

```bash
mkdir -p logs && sbatch hpc/submit_slurm.sbatch
```

## AWS (single EC2 instance + S3)

One-shot instance: boots Amazon Linux 2023, installs Docker, clones this repo,
runs the analysis, syncs outputs to S3, then **terminates itself** (a
`shutdown -h +MAX_HOURS` dead-man switch also guards against hangs).

```bash
export S3_BUCKET=my-bucket          # existing bucket you can write to
export INSTANCE_PROFILE=my-ec2-role # IAM instance profile allowing s3:PutObject
./aws/launch_ec2.sh
```

Default instance `m6i.xlarge` (60 GB gp3). Outputs land in
`s3://$S3_BUCKET/neoantigen-prediction-mhcflurry/<timestamp>/`, with `analysis.log`. Remember to
check cost for your instance type and region; nothing is created until you run
the script, and the instance stops billing when it terminates.

For managed scale-out, `nextflow.config` also has an `awsbatch` profile (fill in your own queue/bucket).

## What was verified

- Nextflow `-profile slurm` / `awsbatch` configs parse (`nextflow config`).
