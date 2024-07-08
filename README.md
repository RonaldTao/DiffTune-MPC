# DiffTune-MPC

This is a DiffTune-MPC toolset for MPC controller auto-tuning using sensitivity propagation. We provide two examples that demonstrate the applications of DiffTune-MPC in matlab to facilitate quick deployment to other applications. 

Details of the DiffTune-MPC can be found in:<br />
[DiffTune-MPC: Closed-Loop Learning for Model Predictive Control]([https://arxiv.org/abs/2209.10021](https://arxiv.org/abs/2312.11384))<br />
[DiffTune: Auto-Tuning through Auto-Differentiation](https://arxiv.org/abs/2209.10021)<br />

If you think this toolset is helpful to your research/application, please cite:<br />
```
@article{tao2024difftune,
  title={DiffTune-MPC: Closed-Loop Learning for Model Predictive Control},
  author={Tao, Ran and Cheng, Sheng and Wang, Xiaofeng and Wang, Shenlong and Hovakimyan, Naira},
  journal={IEEE Robotics and Automation Letters},
  year={2024},
  publisher={IEEE}
}
@article{cheng2022difftune,
  title={DiffTune: Auto-Tuning through Auto-Differentiation},
  author={Cheng, Sheng and Kim, Minkyung and Song, Lin and Wu, Zhuohuan and Wang, Shenlong and Hovakimyan, Naira},
  journal={arXiv preprint arXiv:2209.10021},
  year={2022}
}
```

## Prerequisites

You need to install [https://docs.acados.org/matlab_octave_interface/index.html.](https://docs.acados.org/index.html) on your computer.

Sheng: July 3, tested on Ubuntu 20.04, matlab 2022b, acados github version on July 2, 2024
