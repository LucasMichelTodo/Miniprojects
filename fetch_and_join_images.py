import sys
import os.path
from tqdm import tqdm
from PIL import Image

## Set folders
new = '/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/Plots/'

old = '/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/Plots/'

ratio = 'Ratio/Gene_Level/'
red = 'Red_Signal/Gene_Level/'

outdir = '/mnt/Disc4T/Projects/PhD_Project/Plots/Temp/'
os.makedirs(outdir, exist_ok=True)

## Load genes from a plain text file with one gene ID per line.
# gene_file = '/mnt/Disc4T/Projects/PhD_Project/Microarrays/Gene_lists_for_plots/gens_sense_dif_coverage.txt'
# gids = [line.strip() for line in open(gene_file, 'r+')]

gids = [
    'PF3D7_0102200'
]

## Main Loop
for gid in tqdm(gids, position = 0, leave = True, desc = 'Completed'):

    img_files = [
        new+ratio+gid+'.jpeg',
        old+ratio+gid+'.png',
        new+red+gid+'.jpeg',
        old+red+gid+'.png'
    ]
    
    width, height = 1500, 1500
    images = []
    for im in img_files:
        if os.path.isfile(im):
            images.append(Image.open(im))
        else:
            images.append(Image.new('RGB', (width, height)))
        
    #images = [Image.open(x) for x in img_files]    
    resized_images = [im.resize((width, height)) for im in images]
    
    total_width = 3000
    total_height = 3000
    
    new_im = Image.new('RGB', (total_width, total_height))
    
    new_im.paste(resized_images[0], (0, 0))
    new_im.paste(resized_images[1], (width, 0))
    new_im.paste(resized_images[2], (0, height))
    new_im.paste(resized_images[3], (width, height))

    new_im.save(outdir+gid+'_joined.jpg')
    #print(f'Joining {gid} plots...')


    
