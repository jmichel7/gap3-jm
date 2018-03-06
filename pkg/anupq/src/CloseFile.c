#include <stdlib.h>
#include "pq_defs.h"
#include "constants.h"

/* close file */

void CloseFile (FILE_TYPE file)
{
   if (CLOSE(file) != 0) {
      perror (NULL);
      exit (FAILURE);
   }
}
