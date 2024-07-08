# DiffTune-MPC

This is a DiffTune-MPC toolset for MPC controller auto-tuning using sensitivity propagation. We provide two examples (one linear system and one nonlinear system) that demonstrate the applications of DiffTune-MPC in matlab to facilitate quick deployment to other applications. 

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

You need to install [acados](https://docs.acados.org/index.html) on your computer for matlab interface. Then, clone this repository into the ```/examples``` folder under the installed acados folder.

## Run examples

We offer two examples: a differential wheeled robot and a double integrator system. The system dynamics and constraints are included in ```Differential_Wheeled_Robot.m``` and ```Double_Integrator_System.m```. The code for application of Difftune-MPC are included in ```Differential_Wheeled_Robot_DifftuneMPC.m``` and ```Double_Integrator_System_DifftuneMPC.m```. For the differential wheeled robot, we inlcude a nonlinear system model with linear inequalities on the state and control input. We use acados to solve for the original MPC problem, and use quadprog to solve for the auxilary MPC problems (LMPC-Grad) for the analytical gradients. For the double integrator system, which is linear, we use acaods to solve for both the original MPC problem and the auxilary MPC problems. After running the example, you should be able to see the tracking performance of the closed-loop system using the initial parameters and the learned parameters using Difftune-MPC.

Sheng: July 3, tested on Ubuntu 20.04, matlab 2022b, acados github version on July 2, 2024
