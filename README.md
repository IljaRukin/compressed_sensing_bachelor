# compressed_sensing_bachelor

### license notice

code for my bachelor thesis in compressed sensing (in german) \
the simple solvers IHT (gradient descent with iterative hard threshholding) and OMP (orthogonal matching pursuit) were implemented by me and distributed under GNU General Public License Version 3. \
The other solvers are NESTA (implemented by: J. Bobin, S. Becker; under GNU General Public License 3), FPCAS (implemented by: Zaiwen Wen, Wotao Yin; under GNU GENERAL PUBLIC LICENSE Version 3), spgl1 (implemented by: P. Friedlander, Ewout van den Berg ; under GNU LESSER GENERAL PUBLIC LICENSE Version 2.1). \
The code together with it's license is stored in the folder "solvers".

### Beschreibung

\subsection*{Implementierung}

\subsubsection*{komprimierte Fouriertransformation}
Die Fouriertransformation erzeugt komplexe Koeffizienten, welche für reelle Bilder doppelt besetzt sind. Es wurde die Funktion \textit{compress\_coeff (\textbf{x})} implementiert, welche die Symmetrie der Koeffizienten nutzt um doppelte Koeffizienten zu vernachlässigen und die Koeffizienten in eine reelle Matrix anzuordnen. Die Funktion \textit{decompress\_coeff (\textbf{x})} macht diese Reduktion der Koeffizienten rückgängig.

\subsubsection*{Wavelet-Transformationen}

Es wurden die HAAR und DB4 implementiert. Der Grund dafür war ein Verständnis für die Transformationen zu entwickeln. Außerdem ändern die in Matlab vorhandenen Transformationen die Bildgröße, so dass die vorhandenen Transformationen nicht für CS geeignet waren. Für die Transformation muss ein Signal ($\textbf{f}$) als Zeilenvektor, die Transformation (trafo) als Zeichenkette \glqq haar\grqq oder \glqq db04\grqq und die Ordnung der Transformation (iter) als Zahl übergeben werden. Die Transformation eines eindimensionalen Signals in die Koeffizienten erfolgt mit \textbf{x}=\textit{wave\_1d\_multi (\textbf{f},trafo,iter)} und die Rücktransformation mit \textbf{f}=\textit{iwave\_1d\_multi (\textbf{x},trafo,iter)}. Für die zweidimensionale Transformation wurden zwei verschiedene Schemata umgesetzt:
\begin{itemize}
\item Transformation um eine Ordnung zur Zeit, so dass die Ordnung zeilen- und spaltenweise identisch ist (Funktion: \textit{wave\_2d\_standard (\textbf{X},trafo,iter)})
\item vollständige zeilenweise Transformation und anschließend spaltenweise Transformation. Ordnungen können unterschiedlich sein (Funktion: \textit{wave\_2d\_nonstandard (\textbf{X},trafo,[iter1,iter2])})
\end{itemize}
Dabei muss der Funktion eine Matrix \textbf{X} statt des Vektors übergeben werden.
Für das zweite Schema wird die Transformationsordnung als Vektor mit zwei Zahlen angegeben. Die erste Zahl gibt die Ordnung für die zeilenweise Transformation und die zweite Zahl für die spaltenweise Transformation. Die anderen Argumente sind identisch.
Für alle Transformationen liegen Beispiel-Skripte (test\_1d.m und test\_2d.m) im Ordner Wavelets vor.

\subsubsection*{Vereinfachte Funktion zur Rekonstruktion mit CS}

Für eine einfache Rekonstruktion von Samples wurde eine Funktion erstellt. Die Funktion

\begin{equation*}
\textrm{ [result,fit\_error,comp\_time] = \textit{reconstruct(samples,indices,dim,trafo,solver,sigma)} }
\end{equation*}

nimmt die Samples (samples) mit ihren Positionen (indices) und der Größe des Bildes als Vektor mit der Anzahl an Elementen in einer Zeile und Spalte (dim) entgegen. Des weiteren muss eine Transformation (trafo: cdft, dft2, dct2, haar, db04) und ein Löser (solver: spgl1, nesta, fpcas) als Zeichenkette und das vorhandene Rauschen (sigma) als Zahl angegeben werden. Die Funktion gibt das Ergebnis der Rekonstruktion (result) zusammen mit dem Fehler der Rekonstruktion (fit\_error) und der Rechenzeit (comp\_time) aus. Dazu gibt es ein Beipiel-Skript (final.m). Sollten die implementierten Algorithmen OMP oder IHT verwendet werden, so muss auf die Funktionen direkt zugegriffen werden, weil diese Funktionen mehr Argumente benötigen. Genaueres zu den Parametern ist im Abschnitt zu dem jeweiligen Algorithmen beschrieben. Zu den beiden Algorithmen sind auch Beispiel-Skripte (final\_OMP.m und final\_IHT.m) vorhanden.

Bei all diesen Lösern wird die Matrix als Funktion übergeben. Dies hat den Vorteil, dass die Transformationsmatrix nicht explizit gespeichert werden muss und schneller berechnet werden kann, als eine Matrixmultiplikation.


Der Funktion für die Rücktransformation $\textbf{y} = AT(\textbf{x}) = \textbf{A} \textbf{x}$ werden N x M Koeffizienten $\textbf{x}$ als Vektor übergeben, welche innerhalb der Funktion in eine Matrix umgeformt, transformiert und zurück in einen Vektor angeordnet werden.
Bei diesen Lösern musste auch eine Multiplikation mit der adjungierten Matrix $\textbf{x} = AA(\textbf{y}_n) = \overline{\textbf{A}} \cdot \textbf{y}$ für einen Vektor $\textbf{y}$ formuliert werden.


\subsubsection*{OMP}


\begin{align*}
\textrm{
x = OMP(samples,dim,indices,sparsity,start\_lsqlin\_OptimalityTolerance,final\_lsqlin\_OptimalityTolerance
} \\
\textrm{accuracy,AA,AT)}
\end{align*}

Der Algorithmus OMP nimmt die Samples (samples), Positionen der Samples (indices) und Größe des Bildes (dim) entgegen. Der implementierte Algorithmus wird beendet, sobald die maximale Sparsity (sparsity) oder die geforderte Genauigkeit (accuracy) erreicht wurde. Anschließend werden die endgültigen Koeffizienten mit einer höheren Genauigkeit (final\_lsqlin\_OptimalityTolerance) berechnet. Zusätzlich müssen die Hintransformation und dessen adjungierte (AA und AT) angegeben werden.

\subsubsection*{IHT}

\begin{align*}
\textrm{x = IHT(samples,dim,indices,accuracy,min\_change,maxiter,min\_step,step\_search,} \\
\textrm{start\_sparsity,max\_sparsity,sparsity\_step,AA,AT)}
\end{align*}

Der Algorithmus IHT nimmt die Samples (samples), Positionen der Samples (indices) und Größe des Bildes (dim) entgegen.
In dieser Arbeit implementierte IHT startet mit einer Sparsity von \textit{start\_sparsity} und erhöht diese um \textit{sparsity\_step} bis maximal zur \textit{max\_sparsity}. Zusätzlich müssen die Hin- und Rücktransformationen (AA und AT) angegeben werden.
Der IHT wird abgebrochen, wenn die Fitgenauigkeit $\textrm{\textit{accuracy}} =\sfrac{||\textbf{y}-\textbf{y}^r||_2^2}{||\textbf{y}||_2^2}$ der rekonstruierten Samples $\textbf{y}^r$ überschritten wird. Wenn sich die Fitgenauigkeit um weniger als um \textit{min\_change} ändert, wird die sparsity erhöht und bei erreichen der \textit{max\_sparsity} wird der Algorithmus abgebrochen.\\
Die anfängliche Schrittweite des IHT wurde auf 10 gesetzt. Es wird der Gradient berechnet und ein Schritt in die Richtung größter Abnahme gemacht und der Hard Thresholding Operator auf das Ergebnis angewandt und dann die Fitgenauigkeit berechnet. Die Schrittweite wird um den Faktor \textit{step\_search} bis maximal \textit{min\_step} immer wieder verkleinert solange die Fitgenauigkeit um mindestens \textit{min\_change} zunimmt. Wird \textit{min\_step} unterschritten, so wird die Sparsity erhöht und die Schrittweite erneut gesucht. Ändert sich die Fitgenauigkeit um weniger als \textit{min\_change} oder nimmt sogar ab, so wird der Schritt mit der größten Fitgenauigkeit gemacht.
