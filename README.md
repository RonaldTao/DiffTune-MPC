# DiffTune-MPC

This is a DiffTune-MPC toolset for MPC controller auto-tuning using sensitivity propagation. We provide two examples (one linear system and one nonlinear system) that demonstrate the applications of DiffTune-MPC in matlab to facilitate quick deployment to other applications. 

Details of the DiffTune-MPC can be found in:<br />
[DiffTune-MPC: Closed-Loop Learning for Model Predictive Control](https://arxiv.org/abs/2312.11384)<br />
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

You need to install [acados](https://docs.acados.org/index.html) on your computer for **Matlab + Simulink and Octave Interface**. Then, clone this repository into the ```/examples``` folder under the installed acados folder. Make sure ```env.sh``` is included in the cloned folder.

## Run examples

We offer two examples: a differential wheeled robot and a double integrator system. In both examples, we aim to control the system to track a desired trajectory using MPC. The system dynamics and constraints are included in ```Differential_Wheeled_Robot.m``` and ```Double_Integrator_System.m```. The code for application of Difftune-MPC is included in ```Differential_Wheeled_Robot_DifftuneMPC.m``` and ```Double_Integrator_System_DifftuneMPC.m```. For the differential wheeled robot, we inlcude a nonlinear system model with linear inequalities on the state and control input. We use acados to solve for the original MPC problem, and use quadprog to solve for the auxilary MPC problems (LMPC-Grad) to achieve the analytical gradients. For the linear double integrator system, we use acaods to solve for both the original MPC problem and the auxilary MPC problems. After running the example, you should be able to see two plots showing the tracking performance of the closed-loop system using the initial parameters and the learned parameters using Difftune-MPC.

## Issues/Questions/Suggestions
Feel free to open up an issue if you run into trouble. 

# Authors

**[Ran Tao](https://github.com/RonaldTao)**
**[Sheng Cheng](https://github.com/Sheng-Cheng)**


## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FSheng-Cheng%2FDiffTuneOpenSource&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
