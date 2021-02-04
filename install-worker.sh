# Activate conda
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)";

# Create environment and install dependencies
conda create -y --name mas-worker python=3.8 mysqlclient;
conda activate mas-worker;
conda install -c bioconda -c conda-forge -y glimmer=3.02 blast=2.9.0 trnascan-se=2.0.6 hhsuite;
pip install luigi celery==5.0.2 django==3.1.2 django-simple-history django-debug-toolbar django-crispy-forms djangorestframework biopython==1.77 pandas requests==2.24.0;

# Generate a secret key
python -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())' > "${0%/*}/MAS/settings_files/secret_key.txt"