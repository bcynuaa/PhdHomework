# 目录说明

1. `src`文件夹下存有源代码；
2. `case`文件夹下存有各种程序、格式的算例代码；
3. `image`文件夹下存有各种`gif`图片；
4. `report.pdf`为报告。

以SingleWave方程为例，隐式向后欧拉方法对应的算例文件为`case//SingleWave//implicit_backward.jl`文件，其调用了`src`文件夹下的源码进行计算，而同文件夹下的`animate_ib.jl`为对应的`gif`生成文件。

需要注意的是为了节省上传作业内存，跑完的数据并未上传，但对应的结果在`image`文件夹下的`gif`图片中。如`image//SingleWave//implicit_backward`文件夹下的所有`gif`文件。

同文件夹下的`run_all.sh`脚本为`linux`下的运行脚本，可以运行同目录下所有文件。

另外在`case`文件夹下的`visualize.ipynb`文件为绘图`jupyter-notebook`文件，用于绘制`report.pdf`中的图片。

