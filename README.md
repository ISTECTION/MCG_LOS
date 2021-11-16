# ЛАБОРАТОРНАЯ РАБОТА №3


## Цель работы


Изучить особенности реализации трехшаговых итерационных методов для СЛАУ с разреженными матрицами. Исследовать влияние предобусловливания на сходимость изучаемых методов на нескольких матрицах большой (не менее 10000) размерности.


## Вариант

Сравнить МСГ<sup id="a1">[1](#f1)</sup> и ЛОС<sup id="a2">[2](#f2)</sup> для симметричной матрицы

## О программе

- LOS.hpp
    + [class LOS](#LOS)
- MCG.hpp
    + [class MCG](#MCG)
- generator.hpp
    + [class Generator](#Generator)
    + [class GenerateGilbert](#Generator)
    + [class Generate_Ak](#Generator)
    + [class Generate_diagDomination](#Generator)
- Timer.hpp
    + [class Timer](#Timer)

[#LOS]: fdsfdsfsdfs


<div id="LOS"/>

## LOS (Локально-оптимальная схема)

```c++
none    (); /// метод без предобусловливания
diagonal(); /// метод с диагональным предобуславливанием
hollesky(); /// метод с предобусловливанием Холесского
```


<div id="MCG"/>

## MCG (Метод сопряженных градиентов)

```c++
none    (); /// метод без предобусловливания
diagonal(); /// метод с диагональным предобуславливанием
hollesky(); /// метод с предобусловливанием Холесского
```

<div id="Generator"/>

## Generator

При помощи `class Generator` можно сгенерировать рандомную положительно-определённую матрицу любых размеров

```c++
Generator gen("file/generator/random", 1000, 1E-14, 100000, Chance { 1, 500 });
```
- `"file/generator/random"` - _путь для записи файлов_
- `1000` - _размерность матрицы_
- `1E-14` - _погрешность_
- `100000` - _максимальное количество итераций_
- `Chance { 1, 500 }` - _шанс с которым будут генерироваться элементы в нижнем треугольнике_

При помощи `class GenerateGilbert` можно сгенерировать матрицу Гильберы любух размеров

```c++
GenerateGilbert gen("file/generator/gilbert", 10, 1E-14, 10000);
```
**Параметры аналогичны вышеупомянутым**

При помощи `class Generate_Ak` можно сгенерировать матрицу вида:



<!-- $$
\begin{equation*}
    a_{ii} =
    \begin{cases}
        -\sum\limits_{i \neq j} a_{ij}, i > 1 \\
        -\sum\limits_{i \neq j}a_{ij} + 10^{-k}, i=1
    \end{cases}
\end{equation*}
$$ -->

<img src="https://render.githubusercontent.com/render/math?math=\begin{equation*}a_{ii}=\begin{cases}-\sum_{i \neq j} a_{ij}, i > 1 \\ -\sum_{i \neq j} { a_{ij} } \pm 10^{-k}, i = 1 \end{cases}\end{equation*}"> <br>
**В выражение a_ij + 10^{-k} должен быть** ```+```

<!-- $$A^kx^k=F^k, k = 0, 1, 2, ... $$ -->
<img src="https://render.githubusercontent.com/render/math?math=A^kx^k=F^k, k = 0, 1, 2, ... "> <br>


<!-- $$a_{ii} \in \{ 0, -1, -2, -3, -4 \} $$ -->
<img src="https://render.githubusercontent.com/render/math?math=a_{ii} \in \{ 0, 1, -2, -3, -4 \}"> <br>


```c++
Generate_Ak gen("file/generator/Ak", 10, 1e-14, 100000, 20);
```
- `20` - _k = 0, 1, 2, ... , 19_

**Остальные параметры аналогичны вышеупомянутым**


При помощи `class Generate_diagDomination` можно сгенерировать матрицу с диагональным преобладанием (доминированием):

<!-- $$|a_{ii} | \geq \sum\limits_{j \neq i} |a_{ij}|$$ -->
<img src="https://render.githubusercontent.com/render/math?math=|a_{ii} | \geq \sum_{j \neq i}|a_{ij}|"> <br>



Говорят, что квадратная матрица <!-- $$A_{nn}$$ -->
<img src="https://render.githubusercontent.com/render/math?math=A_{nn}">
обладает свойством диагонального преобладания, если для каждого <!-- $$i=1,\dots,n$$ -->
<img src="https://render.githubusercontent.com/render/math?math=i=1,\dots,n">
причём хотя бы одно из этих неравенств является строгим. Если все неравенства строгие, то говорят, что матрица обладает строгим диагональным преобладанием.


```c++
Generate_diagDomination gen("file/generator/diagonal", 10, 1e-14, 100000);
```
**Параметры аналогичны вышеупомянутым**

<div id="Timer"/>

## Timer

При помощи `class Timer` можно замерять время выполнения методов _следующим образом:_

```c++
Timer::Timer timer;
MCG<double> a("file/generator/1000");
a.solve(Conditional::NONE, false);
timer.setTimeEnd();
std::cout << timer;
```

___

<sub id="f1"> МСГ — Метод сопряженных градиентов [↩](#a1) </sub>

<sub id="f2"> ЛОС — Локально-оптимальная схема [↩](#a2)   </sub>