# compressed_sensing_bachelor

### license notice

code for my bachelor thesis in compressed sensing (thesis is in german) 
the simple solvers IHT (gradient descent with iterative hard threshholding) and OMP (orthogonal matching pursuit) were implemented by me and distributed under GNU General Public License Version 3. 
The other solvers are NESTA (implemented by: J. Bobin, S. Becker; under GNU General Public License 3), FPCAS (implemented by: Zaiwen Wen, Wotao Yin; under GNU GENERAL PUBLIC LICENSE Version 3), spgl1 (implemented by: P. Friedlander, Ewout van den Berg ; under GNU LESSER GENERAL PUBLIC LICENSE Version 2.1). 
The code together with it's license is stored in the folder "solvers".

### Beschreibung

## Implementierung

## komprimierte Fouriertransformation
Die Fouriertransformation erzeugt komplexe Koeffizienten, welche für reelle Bilder doppelt besetzt sind. Es wurde die Funktion ***compress_coeff*** (**x**) implementiert, welche die Symmetrie der Koeffizienten nutzt um doppelte Koeffizienten zu vernachlässigen und die Koeffizienten in eine reelle Matrix anzuordnen. Die Funktion ***decompress_coeff*** (**x**) macht diese Reduktion der Koeffizienten rückgängig.

## Wavelet-Transformationen

Es wurden die HAAR und DB4 implementiert. Der Grund dafür war ein Verständnis für die Transformationen zu entwickeln. Außerdem ändern die in Matlab vorhandenen Transformationen die Bildgröße, so dass die vorhandenen Transformationen nicht für CS geeignet waren. Für die Transformation muss ein Signal (**f**) als Zeilenvektor, die Transformation (trafo) als Zeichenkette "haar" oder "db04" und die Ordnung der Transformation (iter) als Zahl übergeben werden. Die Transformation eines eindimensionalen Signals in die Koeffizienten erfolgt mit **x**=***wave_1d_multi*** (**f**,trafo,iter) und die Rücktransformation mit **f**=***iwave_1d_multi*** (**x**,trafo,iter). Für die zweidimensionale Transformation wurden zwei verschiedene Schemata umgesetzt:
- Transformation um eine Ordnung zur Zeit, so dass die Ordnung zeilen- und spaltenweise identisch ist (Funktion: ***wave_2d_standard*** (**X**,trafo,iter))
- vollständige zeilenweise Transformation und anschließend spaltenweise Transformation. Ordnungen können unterschiedlich sein (Funktion: ***wave_2d_nonstandard*** (**X**,trafo,[iter1,iter2]))
Dabei muss der Funktion eine Matrix **X** statt des Vektors (da 2dimensional) übergeben werden.
Für das zweite Schema wird die Transformationsordnung als Vektor mit zwei Zahlen angegeben. Die erste Zahl gibt die Ordnung für die zeilenweise Transformation und die zweite Zahl für die spaltenweise Transformation. Die anderen Argumente sind identisch.
Für alle Transformationen liegen Beispiel-Skripte (wavelets_test_1d.m und wavelets_test_2d.m).

## Vereinfachte Funktion zur Rekonstruktion mit CS (für NESTA,SPGL1,FPCAS)

Für eine einfache Rekonstruktion von Samples wurde eine Funktion erstellt (wrapper um die ursprünglichen Funktionen). Die Funktion

```
[result,fit_error,comp_time] = ***reconstruct***(samples,indices,dim,trafo,solver,sigma)
```


nimmt die bekannten Datenpunkte, kurz Samples gennant (samples) mit ihren Positionen (indices) und der Größe des Bildes als Vektor mit der Anzahl an Elementen in einer Zeile und Spalte (dim) entgegen. Des weiteren muss eine Transformation (trafo: cdft, dft2, dct2, haar, db04) und ein Löser (solver: spgl1, nesta, fpcas) als Zeichenkette und das vorhandene Rauschen (sigma) als Zahl angegeben werden. Die Funktion gibt das Ergebnis der Rekonstruktion (result) zusammen mit dem Fehler der Rekonstruktion (fit_error) und der Rechenzeit (comp_time) aus. Dazu gibt es ein Beipiel-Skript (CS_final.m). Sollten die implementierten Algorithmen OMP oder IHT verwendet werden, so muss auf die Funktionen direkt zugegriffen werden, weil diese Funktionen mehr Argumente benötigen. Genaueres zu den Parametern ist im Abschnitt zu dem jeweiligen Algorithmen beschrieben. Zu den beiden Algorithmen sind auch Beispiel-Skripte (CS_final_OMP.m und CS_final_IHT.m) vorhanden.

Bei all diesen Lösern wird die Matrix als Funktion übergeben. Dies hat den Vorteil, dass die Transformationsmatrix nicht explizit gespeichert werden muss und schneller berechnet werden kann, als eine Matrixmultiplikation.


Der Funktion für die Rücktransformation **y** = AT(**x**) = **AT** **x** werden N x M Koeffizienten **x** als Vektor übergeben, welche innerhalb der Funktion in eine Matrix umgeformt, transformiert und zurück in einen Vektor angeordnet werden.
Bei diesen Lösern musste auch eine Multiplikation mit der adjungierten Matrix AA **x** = AA(**y**) = ** A^* ** **y** für einen Vektor **y** formuliert werden.


## OMP

```
x = OMP(samples,dim,indices,sparsity,start_lsqlin_OptimalityTolerance,final_lsqlin_OptimalityTolerance,accuracy,AA,AT)
```
Der Algorithmus OMP nimmt die Samples (samples), Positionen der Samples (indices) und Größe des Bildes (dim) entgegen. Der implementierte Algorithmus wird beendet, sobald die maximale Sparsity (sparsity) oder die geforderte Genauigkeit (accuracy) erreicht wurde. Anschließend werden die endgültigen Koeffizienten mit einer höheren Genauigkeit (final_lsqlin_OptimalityTolerance) berechnet. Zusätzlich müssen die Hintransformation und dessen adjungierte (AA und AT) angegeben werden.

## IHT

```
x = IHT(samples,dim,indices,accuracy,min_change,maxiter,min_step,step_search,start_sparsity,max_sparsity,sparsity_step,AA,AT)
```
Der Algorithmus IHT nimmt die Samples (samples), Positionen der Samples (indices) und Größe des Bildes (dim) entgegen.
In dieser Arbeit implementierte IHT startet mit einer Sparsity von *start_sparsity* und erhöht diese um *sparsity_step* bis maximal zur *max_sparsity*. Zusätzlich müssen die Hin- und Rücktransformationen (AA und AT) angegeben werden.
Der IHT wird abgebrochen, wenn die Fitgenauigkeit *accuracy* =||**y**}-**y^r**||_2^2 / ||**y**||_2^2 der rekonstruierten Samples **y^r** überschritten wird. Wenn sich die Fitgenauigkeit um weniger als um *min_change* ändert, wird die sparsity erhöht und bei erreichen der *max_sparsity* wird der Algorithmus abgebrochen.
Die anfängliche Schrittweite des IHT wurde auf 10 gesetzt. Es wird der Gradient berechnet und ein Schritt in die Richtung größter Abnahme gemacht und der Hard Thresholding Operator auf das Ergebnis angewandt und dann die Fitgenauigkeit berechnet. Die Schrittweite wird um den Faktor *step_search* bis maximal *min_step* immer wieder verkleinert solange die Fitgenauigkeit um mindestens *min_change* zunimmt. Wird die minimale Schrittweite *min_step* unterschritten, so wird die Sparsity erhöht und die Schrittweite erneut gesucht. Ändert sich die Fitgenauigkeit um weniger als *min_change* oder nimmt sogar ab, so wird der Schritt mit der größten Fitgenauigkeit gemacht.
