# Simple md support

showcase of **strong** and *em* fonts (or ***both***), links to [Google](google.com)

**Simple** *lists*:

* hi
* hello
* greetings


Also numbered:

1. foo
2. bar
3. stuff

also code support:

```
x = 0
x = 2 + 2
what is x
```

Additionally link to Records, AbstractRecords and Selections. You can also **reference** their fields (for example record key)

link to [[record_root]] and its key: [[record_root#keys#problem]] looks like this

or like this:

- [[record_root#keys#pause_after_run]]
- [[record_root#keys#problem]]
- [[record_SequentialCoupling]]
- [[selection_PartTool]]
- [[record_SoluteTransport_DG_Data#keys#diff_m]]



Using latex in md:

(( \textbf{bold} and \textit{italic} font ))

or like this: (( $x=\frac{1+y}{1+2z^2}$ )) bigger?


(( $$x=\\frac{1+y}{1+2z^2}$$ ))

or

(($$
 \frac{1}{\displaystyle 1+
   \frac{1}{\displaystyle 2+
   \frac{1}{\displaystyle 3+x}}} +
 \frac{1}{1+\frac{1}{2+\frac{1}{3+x}}}
$$))

or

((
\begin{eqnarray*}
 e^x &\approx& 1+x+x^2/2! + \\
   && {}+x^3/3! + x^4/4! + \\
   && + x^5/5!
\end{eqnarray*}
))

ca
asc