# A simple shell script to test the PyFiNeR fitting

echo "The output should be:"
echo "4779 2  0.5911760 15.166 15.164 16.472 16.468 16.469 16.465 15.579 15.577 15.574 15.572  0.698841  0.066634 -0.019603 -0.059386 0.0323 0.0544 0.000504  108  108"
echo " , while the 4779_2.pdf should look the same as the included 4779_2_comp.pdf"

python ../pyfiner.py "4779 2" 0.591176 4779.k2 4779.j2 4779.h2 4779_2.pdf
