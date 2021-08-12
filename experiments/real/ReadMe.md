**This is a very very basic Readme. Will be updated soon!**



**Step 1.** clone the full repo & submodules

**Step 2.** go to  `./experiments/real` to launch the job

**Step 3.** launch these 2 commands:

`module load anaconda/python3`

`./pyimaging.py params.csv`

**INPUT** in `.csv`

1. imagename

2. run id

3. subcube id 

4. algo

5. gam 

6. gambar

7. rw

8. nReweights

**PS:** data files and G matrices from calib pre-processing step are in ` /lustre/home/shared/sc004/FACETED_HYPERSARA_EXPERIMENTS` . Do not create copies as these are read only and the code is already directed to them. 

**PPS:** new matlab functions are in `./experiments/real/real_data`

