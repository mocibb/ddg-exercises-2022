# ddg-exercises

## 更新
在wiki中陆续加入作业速通的笔记 [作业速通](https://github.com/mocibb/ddg-exercises-2022/wiki/%E4%BD%9C%E4%B8%9A%E9%80%9F%E9%80%9A)

## 正文

这个项目是Keenan Crane教授在CMU开设的离散微分几何的作业的答案。2022年的课程的主页在[Discrete Differential Geometry](https://brickisland.net/DDGSpring2022/) (15-458/858). 

这门课程主要介绍了几何处理相关的数学工具，用这些工具完成有意思的实验。

通过跑通和完成这些实验可以在边学边玩中掌握抽象晦涩的数学概念。

这门课程有一定难度，初学者不理解笔记中的概念或者没有思路完成作业时可以先下载代码跑通实验。

一边实验一边学习，这样能提高学习效率，通过理解作业也让学习更有针对性，也是我把答案放到网上的初衷。

如果遇到不明白的地方，可以不用犹豫添加我的微信 mocibb， 我已经给很多人讲解过DDG中的数学和算法。

在这门课程中我们会学到，单纯复形，网格的半边表示，微分形式和外微分，微分几何，几何处理中的重要工具Laplace-Beltrami算子，曲面的共形参数化方法，向量场的设计等。

这门课程的实验既包括经典的几何处理的算法（求解Poisson方程，计算曲率流等），也有最近几年的图形学顶刊的算法（计算测地线的热核算法，用来设计向量场的平直联络算法）。

通过这门课既可以学习到有趣的几何处理算法，也能通过实验摸到现代数学的门槛，包括代数拓扑，黎曼几何，Hodge定理等。

当然这里面最核心的是离散微分形式，微分形式的语言已经深深地影响到现代的数学和物理，离散微分几何是学习微分形式的最佳路径。

这门课是学习DDG最好的资源，Keenan Crane教授从2016年就开始把授课内容放到网上方便大家学习，这里对Keenan Crane教授表示由衷的感谢。

除了这门课还有一些优秀的学习资源可以参考，感谢各位教授的无私奉献。

- 这里有SIGGRAPH中历次介绍DDG的笔记，可以从最新的看起 http://ddg.cs.columbia.edu/
- UCSD Albert Chern教授的[CSE 270](https://cseweb.ucsd.edu/~alchern/teaching/cse270_wi24/)，Albert Chern教授的Houdini真是出神入化，Slide也非常漂亮
- MIT Justin Solomon教授的[6.8410](https://groups.csail.mit.edu/gdpgroup/68410_spring_2023.html)，Justin Solomon教授的课程有视频，数学推导炉火纯青

## 在线体验版本
  Keenan Crane教授在网上公开了这门课作业的在线体验版本 [geometry-processing-js](https://geometrycollective.github.io/geometry-processing-js/)

   <img src="https://github.com/mocibb/ddg-exercises-2022/assets/18642/f0b13a39-61fa-4285-9879-8f56456d3acd" width="640">
  包括平均曲率流，离散曲率和法向，计算测地线的热核算法，共形参数化方法，向量场的Hodge分解，求解曲面上的Poisson方程，离散微分形式，最后是设计向量场的平直联络算法。
   


## 如何跑通第一个DDG例子。
1. 克隆这个项目，注意这个项目里面有子模组，需要使用--recursive选项。projects里面就是需要完成的作业，我们用geodesic-distance作为例子。
```
git clone --recursive https://github.com/mocibb/ddg-exercises-2022
cd ddg-exercises-2022/projects/geodesic-distance
```
2. 编译
```
mkdir build
cd build
cmake ..
make -j4
```
3. 运行
```
bin/main
```
<img src="https://github.com/mocibb/ddg-exercises-2022/assets/18642/5ec1c0cc-4192-40b6-8c3e-a000aaee2167" width="640">

<img src="https://github.com/mocibb/ddg-exercises-2022/assets/18642/186035e8-86eb-4f6c-af56-d78fe73449b9" width="640">


## License

[MIT](https://opensource.org/licenses/MIT)
