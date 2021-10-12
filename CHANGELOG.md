## VERSION 1.1 - 4/20/21
- Added configuration option to set rabbitmq password.
- Added configuration options to set work log and pid file locations.
- Adding configuration to set conda executable location so mas-worker can be called without sourcing user's bash_rc.
- Luigi configuration file is now linked to mas container as volume, allowing configuration changes to be reflected immediately in the container.
- Added mas-worker.sh script to replace run-worker.sh. This scrip can be used by a daemon.
- Added mas-worker.service to serve as a template for daemonization of the mas worker on machines with systemd.
- Updated README to fix issues with database backup commands.
- Added instructions for daemonizing mas worker to README.
- MAS container now runs Luigi daemon as unprivileged user.
- Fixing bug where submitting searches for all proteins in a genome would not work if on non-CDS feature.
- Fixing hhsearch/PDB evalue formatting.
- Fixing incorrect links to PDB when chain id is 2 letters.
- Fixing result table highlighting.

## VERSION 1.2 - 5/12/21
- Fixed missing secret key bug in the MAS container dockerfile.
- Set first residue to M in protein sequences with non-standard start codons.

## VERSION 1.3 - 10/08/21 (updated 10/12/21)
- Fixed vulnerability: The broker URL for mas-worker is no longer set within the celery multi command. Instead it is set in the settings file.
- Updated httpd to latest version
- Modified the way sacct results are parsed to account for 'extra information' in the job status.
- When uploading a phage with DTRs, search for identical repeat region annotation in database and use if found. This is to prevent a duplicate entry error upon genome upload for phages with the same DTR as a previously uploaded phage.
- Added a button to the result viewer which copies the protein sequence to the user's clipboard.
- Fixed bug causing alignment text to be associated with the incorrect table row in certain instances.
- The 'too many features' warning is now only displayed on the genome navigator.
- MAS will now detect if a worker is running and if not, it will warn the user with text in the bootstrap navbar and disable search and bacterial genome submit buttons.
- Fixed 'annotation' misspelling in genome visualization tooltip
- Fixing incorrect capitalization of 'always' in mas-worker.service which was preventing the service from restarting
- Refactored code, renamed and moved some files, and removed unused code.
- Fixed bug in initial release on 10/08/21: Javascript error caused results not to show if a search has not been previously ran