# 框架

[toc]

## case模块

1. T_inf = 1.
2. C_inf = 1.
3. rho_inf = 1.4
4. P_inf = 1.
5. gamma = 1.4
6. R = 1. / 1.4;
7. Cv = 2.5 / 1/4
8. Ma_inf: 用户输入
9. theta_inf: 用户输入
10. k2 = 0.8（二阶人工粘性）
11. k4 = 0.003（四阶人工粘性）
12. CFL = 2.0
13. error = 1e-5
14. STEP = 10000
15. data_path: 存数据的路径
16. image_path: 存图片的路径

```julia
function expand([rho rho*U rho*V rho*E])
    return [rho U V E P H C]
end
```

## grid模块

读取网格文件

## Flux模块

用于求解每一步计算

## RK4模块

龙格库塔时间推进

## post模块

后处理模块

## EulerEq模块

主程序

需要用户输入来流Ma和迎角alpha
