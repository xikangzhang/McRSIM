# Multi-Camera Robust Shape Interaction Matrix (McRSIM)

This code is implementation for the following paper:

Xikang Zhang, Bengisu Osbay, Mario Sznaier, Octavia Camps, Dynamics Enhanced Multi-Camera Motion Segmentation from Unsynchronized Videos, ICCV 2017.

## Prerequisites

To get experimental results of the paper, you need to download the following data set:

- Hopkins 155 data set: http://www.vision.jhu.edu/data/hopkins155/
- Hopkins Additional Sequences with Missing data: http://vision.jhu.edu/data/
- RSL 12 data set: http://robustsystems.coe.neu.edu/?q=content/publications

## Testing

The demo code are in matlab folder. You may need to set up the correct path before running the code.

- testBenchmark.m: simulated Multi-camera Motion Segmentation on modified Hopkins 155 data set
- testBenchmark2.m: Multi-camera Motion Segmentation on RSL 12 data set
- testMissing.m: Single-camera Motion Segmentation on Hopkins 12 Real Motion Sequences with Incomplete data.
- testMissingMC.m: Multi-camera Motion Segmentation on Hopkins 12 Real Motion Sequences with Incomplete data, half trajectories rotated 45 degree as Camera 2.
- testGross.m: Single-camera Motion Segmentation on Hopkins 155 data set with gross contamination
- testGrossMC.m: Multi-camera Motion Segmentation on Hopkins 155 data set with gross contamination, half trajectories rotated 45 degree as Camera 2.
- testSingleCam.m: Single-camera Motion Segmentation on Hopkins 155 data set

## Citation

Please cite the paper if it helps your research:

    @inproceedings{zhang2017mcms,
        author = {Xikang Zhang, Bengisu Osbay, Mario Sznaier, Octavia Camps},,
        booktitle = {ICCV},
        title = {Dynamics Enhanced Multi-Camera Motion Segmentation from Unsynchronized Videos},
        year = {2017}
    }
