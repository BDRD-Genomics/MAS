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