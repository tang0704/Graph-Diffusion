# MMSkeleton

## Introduction

MMSkeleton is an open source toolbox for skeleton-based human understanding.
It is a part of the [open-mmlab](https://github.com/open-mmlab) project in the charge of [Multimedia Laboratory, CUHK](http://mmlab.ie.cuhk.edu.hk/).
MMSkeleton is developed on our research project [ST-GCN](https://github.com/yysijie/st-gcn/blob/master/OLD_README.md).

<p align="center">
    <img src="demo/recognition/demo_video.gif", width="700">
</p>

## Updates
[2019-08-29] MMSkeleton v0.5 is released.

## Features

- **High extensibility**

    MMSkeleton provides a flexible framework for organizing codes and projects systematically, with the ability to extend to various tasks and scale up to complex deep models.

- **Multiple tasks**

    MMSkeleton addresses to multiple tasks in human understanding, including but not limited to:
    - skeleton-based action recognition: [[ST-GCN]](./doc/START_RECOGNITION.md)
    - skeleton-based action generation
    - 2D/3D pose estimation
    - pose tracking

## Getting Started
Please see [INSTALL.md](./doc/INSTALL.md) and [START_RECOGNITION.md](./doc/START_RECOGNITION.md) for the basic usage of MMSkeleton.

## License
The project is release under the [Apache 2.0 license](./LICENSE).

## Contributing
We appreciate all contributions to improve MMSkeleton.
Please refer to [CONTRIBUTING.md](./doc/CONTRIBUTING.md) for the contributing guideline.


## Citation
Please cite the following paper if you use this repository in your reseach.
```
@inproceedings{stgcn2018aaai,
  title     = {Spatial Temporal Graph Convolutional Networks for Skeleton-Based Action Recognition},
  author    = {Sijie Yan and Yuanjun Xiong and Dahua Lin},
  booktitle = {AAAI},
  year      = {2018},
}
```

## Contact
For any question, feel free to contact
```
Sijie Yan     : ys016@ie.cuhk.edu.hk
Yuanjun Xiong : bitxiong@gmail.com
```
