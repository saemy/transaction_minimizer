Transaction minimizer
================

When being on holydays it might happen that people are paying stuff that is also
for others. At the end these debts have to be settled up. This tool is doing so
by moreover minimizing the overall number of transactions that are needed.

Usage
---------

The tool reads the input directly from stdin:

```shell
./transaction_minimizer <<EOF
Leo Alice 5
Leo Bob 9
Leo Theo 4
Leo Adam 7
Lucy Alice 1
Lucy Bob 3
Alice Bob 2
Alice Theo 6
Theo Adam 8
Rachel Lucy 10
Rachel Bob 12
Rachel Theo 11

Ivy Timmy 100
Timmy Ivy 40
EOF
```

Every line consists of three columns: The first lists the person who owes the
one in the second the amount in the third. You can add as many spaces or tabs
between the columns as you like.

It then prints the transactions to stdout:

```shell
Theo   -> Adam   : 15.00
Ivy    -> Timmy  : 60.00
Alice  -> Bob    : 26.00
Alice  -> Theo   : 28.00
Leo    -> Alice  : 25.00
Lucy   -> Alice  : 27.00
Rachel -> Lucy   : 33.00
```
this equals the following payment graph
```shell
Ivy --60--> Timmy

Rachel --33--> Lucy --27-->       /--28--> Theo --15--> Adam
                            Alice
                Leo --25-->       \--26--> Bob
```

Note, since we reduce as much transactions as possible, it is likely that a
 person only has to pay one big amount to another one. This guy then forwards it
to the next person until the destination is reached. This requires the people to
a) trust each other that the payments are really forwarded and b) payments may
need to be timed s.t. a person in the middle of the payment graph (Alice) has
enough money to do her transactions which contain all forwarded payments.

How it works
------------
The algorithm works in different phases:

1. All cyclic payments are removed. This is done by subtracting the smallest
   amount along a circle from all others and then removing its edge:

  Assume w.l.o.g. ![equation](http://www.sciweavers.org/tex2img.php?eq=p_n%20%3D%20%5Cmin_%7Bi%3D%5C%7B1..n%5C%7D%7D%28p_i%29&bc=Transparent&fc=Black&im=png&fs=12&ff=arev&edit=0)
  then ![equation](http://www.sciweavers.org/tex2img.php?eq=A_1%20%5Coverset%7Bp_1%7D%7B%5Clongrightarrow%7D%20A_2%20%5Coverset%7Bp_2%7D%7B%5Clongrightarrow%7D%20...%20%5Coverset%7Bp_%7Bn-1%7D%7D%7B%5Clongrightarrow%7D%20A_n%20%5Coverset%7Bp_n%7D%7B%5Clongrightarrow%7D%20A_1%20%5CLongrightarrow%20A_1%20%5Coverset%7Bp_1%20-%20p_n%7D%7B%5Clongrightarrow%7D%20A_2%20%5Coverset%7Bp_2%20-%20p_n%7D%7B%5Clongrightarrow%7D%20...%20%5Coverset%7Bp_%7Bn-1%7D%20-%20p_n%7D%7B%5Clongrightarrow%7D%20A_n&bc=Transparent&fc=Black&im=png&fs=12&ff=arev&edit=0). 

2. The payment span is computed. Since all circles are removed by now we operate
   on a directed acyclic payment graph (DAG). However, it can still occur that
   it contains unnecessary transactions:
   ```
   A --------10------> B       ===>       A --20--> C --20--> B
    \--10--> C --10-->
   ```
   The span is found by doing a special form of a topological sort. It starts at
   the nodes that only receive money (the leaves in the DAG) and then operates
   in rounds: In round *i* all nodes that are reachable in *i* hops from the
   leaves (using reversed edges) are inspected. Nodes that already are in the
   span are ignored. If an observed node only has edges to nodes that already
   are in the span, it is added. Moreover, the edges that link to nodes that
   were added in round *i-1* are added as well (There allways is at least one
   such edge since otherwise the node would have been added in a previous round).
   Eventually, all nodes are included in the span and this phase is finished.

3. The payments along the edges in the span are computed. This is done by also
   building a *routing table* in the previous step. That table tells for a given
   node which outgoing edge that has to be used to reach some destination node.
   Payment amount calculation is started at the nodes that only have to pay (the
   roots of the span). They are marked with the destination node and sent along
   the according edge. A breath-first search propagates the payments through the
   payment span.

TODOs
-----
- Check if the outgoing- minus the incoming payments equal at each node in the
  original graph and the payment span.
- If multiple spans are available with the same number of transactions, choose
  the one with the smallest (maximum) transaction amount.
- Add a flag if payment tunneling should be used or if only circlic payments
  should be removed.
