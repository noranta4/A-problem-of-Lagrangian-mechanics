# A problem of Lagrangian mechanics
The motion of massive point with many constraints. Application of Runge-Kutta methods to find approximate solutions of equations of motion.

University project • 2014 - Laboratorio di Fisica Computazionale - BSc in Physics, II year

The full relation can be found in `A_problem_of_Lagrangian_mechanics.pdf`. The code is in `A_problem_of_Lagrangian_mechanics.c`
## Statement of the problem (italian)

![problem](https://raw.githubusercontent.com/noranta4/A-problem-of-Lagrangian-mechanics/master/img/problem.PNG)

In un piano verticale è posto un semidisco
**d** omogeneo pesante di massa md. Tale
semidisco è ottenuto da un disco di centro
O e raggio r tagliato lungo un diametro AB,
è indicato con G il baricentro del
semidisco. Il punto O è fisso ed il
semidisco può ruotare senza attrito
attorno ad O. Solidalmente ad AB è posta
una guida liscia di massa trascurabile.
Lungo tale guida si muove senza attrito un
punto pesante **p** di massa mp, soggetto
oltre alla forza peso e alla reazione
vincolare alla forza elastica F= -k(Op) con
k >0.

## Running the program (italian)

Il programma chiede di scegliere le costanti del sistema attraverso il canale standard di input. È
possibile scegliere due configurazioni preimpostate con λ>1 o λ<=1, sceglierle a piacere o sceglierle in
modo casuale fra 0 e 100.
In seguito è richiesto di specificare le condizioni iniziali desiderate, tale scelta dipende dalla
precedente e sarà possibile nuovamente sceglierle a piacere, sceglierle in modo casuale fra 0 e 1000
con θ che varia invece fra 0 e 2π, o scegliere due particolari configurazioni con oscillazioni in fase o in
controfase.
Viene chiesto infine di specificare il dt e il tempo totale di integrazione. Il programma stampa su file
sempre lo stesso numero di punti (~10000) affinché il file di output sia di circa 1,5MB, in questo modo
viene preservata la memoria e contenuto notevolmente il tempo di elaborazione (la scrittura su file è
uno dei processi più lenti). Più il tempo di integrazione è lungo più è grande l’intervallo di tempo fra
due punti consecutivi sul file.
I file di output sono nominati `anCOST%CI%.txt` dove al posto di % ci sono le lettere corrispondenti
alla particolare scelta di costanti e condizioni iniziali.

## Examples of **p** trajectories in cartesian coordinates and a plot of the potential energy **U** vs **x** and **theta**.
![flower](https://raw.githubusercontent.com/noranta4/A-problem-of-Lagrangian-mechanics/master/img/flower.PNG)

![net](https://raw.githubusercontent.com/noranta4/A-problem-of-Lagrangian-mechanics/master/img/net.PNG)

![energy](https://raw.githubusercontent.com/noranta4/A-problem-of-Lagrangian-mechanics/master/img/energy.PNG)

