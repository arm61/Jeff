# Jeff - Julian Reflectometry

Jeff is designed to be a **simple** reflectometry analysis package. By **simple** I mean that Jeff does not have any user-friendly bells and whistles, but can access a very powerful backend.

In Jeff, you will **not** find a graphical user-interface for creating models or sliders to modify to your model. The aim of Jeff is to leverage the high performance and broad set of [data science relevant tools](https://turing.ml/dev/) in the Julia language to produce a very fast and powerful reflectometry analysis. This may include auto-differentation libraries and GPU optimisation (although the latter [may not be useful](https://gist.github.com/andyfaff/6f4c8a17fb8f9af804f395b932e6cac0)).

The ambition is that Jeff will be efficient enough to run **during** a given experiment (as the data is collected), and therefore enable the user to access information such as inverse uncertainties and potential Bayesian model evidence in near real time. This will use the powerful samplers of the [Turing.jl](https://turing.ml/dev/) package as well as the live visualisation of [TuringCallbacks](https://turinglang.github.io/TuringCallbacks.jl/dev/). 

In order to achieve this, the following must first be achieved: 

- Complete optimisation of the reflectometry calculation in Julia 
- Assessment of the code to run with the Turing samplers
- Parameter constraining
- Multiple contrast co-refinement
- Porting the Julian reflectometry code to run on GPUs 

If using this, it is possible to obtain some description of the parameter uncertainties in a reasonable time-frame (seconds to minutes timescales), then it may be possible to use this package during a given experiment to enable **experiment-driven science**. 



