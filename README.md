# DiffTune-MPC

This is a DiffTune-MPC toolset for auto-tuning model predictive controllers (MPC) using sensitivity propagation. We provide two examples (one linear system and one nonlinear system) that demonstrate the applications of DiffTune-MPC in MATLAB to facilitate quick deployment to other applications. 

Details of the DiffTune-MPC can be found in:<br />
[DiffTune-MPC: Closed-Loop Learning for Model Predictive Control](https://arxiv.org/abs/2312.11384)<br />
[DiffTune: Auto-Tuning through Auto-Differentiation](https://arxiv.org/abs/2209.10021)<br />

If you think this toolset is helpful to your research/application, please cite:<br />
```
@article{tao2024difftuneMPC,
  title={{DiffTune-MPC}: Closed-Loop Learning for Model Predictive Control},
  author={Tao, Ran and Cheng, Sheng and Wang, Xiaofeng and Wang, Shenlong and Hovakimyan, Naira},
  journal={IEEE Robotics and Automation Letters},
  year={2024},
  publisher={IEEE}
}
@article{cheng2024difftune,
  title={{DiffTune}: Auto-Tuning through Auto-Differentiation},
  author={Cheng, Sheng and Kim, Minkyung and Song, Lin and Yang, Chengyu and Jin, Yiquan and Wang, Shenlong and Hovakimyan, Naira},
  journal={IEEE Transactions on Robotics},
  year={2024}
}
```

## Prerequisites

You need to install [acados](https://docs.acados.org/index.html) on your computer and also install the [Matlab + Simulink and Octave Interface](https://docs.acados.org/matlab_octave_interface/index.html). Then, clone this repository into the ```ACADOS_DIR/examples/acados_matlab_octave``` folder, where ACADOS_DIR is the directory of acados. Make sure to source ```env.sh``` under ACADOS_DIR/examples/acados_matlab_octave/DiffTune-MPC before starting MATLAB. The code has been tested in Ubuntu 20.04.

## Run examples

We offer two examples: a differential wheeled robot and a double integrator system. In both examples, we aim to control the system to track a desired trajectory using MPC. 

For the double integrator system, the linear dynamics and constraints are included in ```Double_Integrator_System.m```. The code for the application of Difftune-MPC is included in ```Double_Integrator_System_DifftuneMPC.m```. We use acaods to solve both the original MPC problem and the auxiliary MPC problems. After running the system for 20 iterations, you should be able to see a figure showing the tracking performance of the closed-loop system using the initial parameters and the learned parameters obtained by Difftune-MPC.

```
$ cd ACADOS_DIR/examples/acados_matlab_octave/DiffTune-MPC 
$ source env.sh
$ matlab
Inside MATLAB, run 
run('Double_Integrator_System_DifftuneMPC.m');
```

For the differential wheeled robot, the nonlinear dynamics and constraints are included in ```Differential_Wheeled_Robot.m```. The code for the application of Difftune-MPC is included in ```Differential_Wheeled_Robot_DifftuneMPC.m```. We include a nonlinear system model with linear inequalities on the state and control input. We use acados to solve the original MPC problem, and use quadprog to solve the auxilary MPC problems (LMPC-Grad) to compute the analytical gradients. After running the system for 20 iterations, you should be able to see a figure showing the tracking performance of the closed-loop system using the initial parameters and the learned parameters obtained by Difftune-MPC.

```
$ cd ACADOS_DIR/examples/acados_matlab_octave/DiffTune-MPC 
$ source env.sh
$ matlab
Inside MATLAB, run 
run('Differential_Wheeled_Robot_DifftuneMPC.m');
```

## Issues/Questions/Suggestions
Feel free to open up an issue if you run into trouble. 

# Authors

**[Ran Tao](https://github.com/RonaldTao)**
**[Sheng Cheng](https://github.com/Sheng-Cheng)**


## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details

[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FRonaldTao%2FDiffTune-MPC&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
