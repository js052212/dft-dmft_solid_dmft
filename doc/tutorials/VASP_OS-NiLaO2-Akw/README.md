这不是官网的例子，是个人在使用过程中的使用例子，
主要目的是使用VASP+solid 绘制谱函数， 因为官网没有相关的后处理的流程，所以本例子中的核心代码都是自己写的，不是官网的，使用者可根据自己使用情况进行一定的修改，但是核心逻辑是没有问题的。 
Alternatively, currently you can do spectral function plots along a kpath only using the w90 interface when having a _hr.dat file one disk. 
Then you can do it via post processing routines (pcb module) from solid_dmft along any path. 
以上是大致流程和文件，刚详细信息请看reference.
