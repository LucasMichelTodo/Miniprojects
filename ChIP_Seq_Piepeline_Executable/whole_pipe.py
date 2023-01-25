import subprocess as sp

sp.run(['python3', 'preprocess_and_align.py'])
sp.run(['python3', 'remove_duplicates.py'])
sp.run(['python3', 'make_tracks.py'])
