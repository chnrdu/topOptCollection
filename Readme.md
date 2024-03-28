# topOptCollection Readme

The collections for Structural Topology Optimization from literatures

The copyright of original programs belong to original authors or publishers. The copyright information or original published paper title is included in the comments of the programs.

The comments are copyrighted by original comments' authors; some comments are updated by Run Du and his team.

从文献中获得的结构拓扑优化算法代码合集。

原始程序的版权归原始作者或出版社所有。原始程序注释中包括原始作者或出版社的版权信息，原始论文标题也包含在注释中。

原始注释的版权归原始注释作者所有；一些注释被杜润及他的团队更新。

Examples:
以下为代码的示例：

1. top99.m
   example:   

    ```
    top99(100,50,0.3,3,3)
    ```

2. top88.m
  example:   

    ```
    top88(100,50,0.3,3,3,1)
    ```

3. multitop.m
  example:   

    ```
    pkg load image % if use octave
    [nx,ny,tol_out,tol_f,iter_max_in,iter_max_out,p,q,e,v,rf] = multitop_setparameters()
    multitop(nx,ny,tol_out,tol_f,iter_max_in,iter_max_out,p,q,e,v,rf)
    ```

3. multitop_h.m
  example: 
  ```
    pkg load image % if use octave
    [nx,ny,tol_out,tol_f,iter_max_in,iter_max_out,p,q,e,v,rf] = multitop_setparameters()
    multitop_h(nx,ny,tol_out,tol_f,iter_max_in,iter_max_out,p,q,e,v,rf)
  ```

Code test platform 代码测试平台

|Code function | MATLAB | octave|
|:---:|:---:|:---:|
| top99 |  ✔ | ✔ |
| top88 | ✔ | ✔ |
| multitop | ✔ | ✔ |
| multitop_h | ✔ | ✔ |