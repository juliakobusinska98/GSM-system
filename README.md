# GSM-system
C++ program implementation of the convolutional code system with the efficiency of 1/6 for TCH / F2.4 data in the GSM system - Viterbi hard-decision decoder
The convolutional code is a redundant code. In order to protect the data against errors, the information sent is accompanied by additional redundant bits, which are used in the receiver to determine the correctness of the received data. The idea of convolutional encoding is to transform the input k-bit information string into n-bit
output string.
Mentioned convolutional generates 6 output bits in response to 1 input bit and is using a 4-block register.
link to the encoder documentation from the task:
https://www.etsi.org/deliver/etsi_en/300900_300999/300909/08.0501_60/en_300909v080501p.pdf
