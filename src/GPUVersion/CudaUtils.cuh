#include <stdio.h>

// La macro "avvolge" la funzione e cattura il NOME DEL FILE e la RIGA
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

// La funzione vera e propria che analizza l'errore
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      // Se c'è un errore, stampa la stringa esatta dell'errore, il file e la riga
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      // E fa crashare il programma di proposito invece di proseguire coi dati corrotti
      if (abort) exit(code);
   }
}