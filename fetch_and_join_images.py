import sys
import os.path
from tqdm import tqdm
from PIL import Image

## Set folders
new = '/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/Plots/'

old = '/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/Plots/'

ratio = 'Ratio/Gene_Level/'
red = 'Red_Signal/Gene_Level/'

snapshots = '/mnt/Disc4T/Projects/Miniprojects/Whole_Genome_IGV_Snapshots/All_gene_snapshots/'

outdir = '/mnt/Disc4T/Projects/PhD_Project/Whole_Genome_Plots/'
os.makedirs(outdir, exist_ok=True)

## Load genes from a plain text file with one gene ID per line.
gene_file = outdir+'all_gene_ids.txt'
gids = [line.strip() for line in open(gene_file, 'r+')]

# gids = [
#     'PF3D7_0102200'
# ]

## Main Loop
for gid in tqdm(gids, position = 0, leave = True, desc = 'Completed'):

    img_files = [
        new+ratio+gid+'.jpeg',
        old+ratio+gid+'.png',
        new+red+gid+'.jpeg',
        old+red+gid+'.png',
        snapshots+'igv_snapshot_'+gid+'.png'
    ]
    
    width, height = 1500, 1500
    images = []
    for im in img_files[0:4]:
        if os.path.isfile(im):
            images.append(Image.open(im))
        else:
            images.append(Image.new('RGB', (width, height)))

    if os.path.isfile(img_files[4]):
        images.append(Image.open(img_files[4]))
    else:
        images.append(Image.new('RGB', (width*2, height*2)))

        
        
    #images = [Image.open(x) for x in img_files]    
    resized_images = [im.resize((width, height)) for im in images[0:4]] + [images[4].resize((width*2, height*2))]
    
    total_width = 6000
    total_height = 3000
    
    new_im = Image.new('RGB', (total_width, total_height))
    
    new_im.paste(resized_images[0], (0, 0))
    new_im.paste(resized_images[1], (width, 0))
    new_im.paste(resized_images[2], (0, height))
    new_im.paste(resized_images[3], (width, height))
    new_im.paste(resized_images[4], (width*2, 0))


    new_im.save(outdir+gid+'_joined.jpg')
    #print(f'Joining {gid} plots...')


    
