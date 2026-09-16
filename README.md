# CrisprTilingTool
CrisprTilingTool is a tool to to detect all possible guides in a gene and report these 




docker build --platform linux/amd64 -t shiny_crispr .
docker run --platform linux/amd64 -d -p 3838:3838 shiny_crispr

## For Google Cloud
docker tag shiny_crispr europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/shiny_crispr:1.0
docker push europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/shiny_crispr:1.0

## Google Cloud Run

gcloud run deploy shiny-crispr \
    --image europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/shiny_crispr \
    --network 'projects/pmc-vpc-res-private-20gx/global/networks/shared-vpc-res-priv-dev' \
    --subnet 'projects/pmc-vpc-res-private-20gx/regions/europe-west4/subnetworks/subnet-res-priv-dev' --service-account='sa-nextflow-runner@pmc-gcp-box-d-pip-development.iam.gserviceaccount.com' \
    --vpc-egress=private-ranges-only

gcloud run services proxy shiny-crispr --project pmc-gcp-box-d-pip-development --region europe-west4