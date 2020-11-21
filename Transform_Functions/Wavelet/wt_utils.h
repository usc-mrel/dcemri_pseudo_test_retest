
#ifndef WT_UTILS_H
#define WT_UTILS_H

#define MACRO_CONCAT(a,b) a##b

#define CONVDWN_FUNC_NAME(ARG_NAME) MACRO_CONCAT(convdwn_,ARG_NAME)

#define TYPE_DIM1 
   #define TYPE_REAL
      #define TYPE_DOUBLE
         #define TYPE_NAME 1rd
         #include "fwt_utils.c"
      #undef TYPE_DOUBLE
         #define TYPE_NAME 1rf
         #include "fwt_utils.c"        
   #undef TYPE_REAL     
      #define TYPE_DOUBLE
         #define TYPE_NAME 1cd
         #include "fwt_utils.c"
      #undef TYPE_DOUBLE
         #define TYPE_NAME 1cf
         #include "fwt_utils.c"
#undef TYPE_DIM1
   #define TYPE_REAL
      #define TYPE_DOUBLE
         #define TYPE_NAME 2rd
         #include "fwt_utils.c"
      #undef TYPE_DOUBLE
         #define TYPE_NAME 2rf
         #include "fwt_utils.c"        
   #undef TYPE_REAL     
      #define TYPE_DOUBLE
         #define TYPE_NAME 2cd
         #include "fwt_utils.c"
      #undef TYPE_DOUBLE
         #define TYPE_NAME 2cf
         #include "fwt_utils.c"        
        
        


#define UPCONV_FUNC_NAME(ARG_NAME) MACRO_CONCAT(upconv_,ARG_NAME)

#define TYPE_DIM1 
   #define TYPE_REAL
      #define TYPE_DOUBLE
         #define TYPE_NAME 1rd
         #include "iwt_utils.c"
      #undef TYPE_DOUBLE
         #define TYPE_NAME 1rf
         #include "iwt_utils.c"        
   #undef TYPE_REAL     
      #define TYPE_DOUBLE
         #define TYPE_NAME 1cd
         #include "iwt_utils.c"
      #undef TYPE_DOUBLE
         #define TYPE_NAME 1cf
         #include "iwt_utils.c"
#undef TYPE_DIM1
   #define TYPE_REAL
      #define TYPE_DOUBLE
         #define TYPE_NAME 2rd
         #include "iwt_utils.c"
      #undef TYPE_DOUBLE
         #define TYPE_NAME 2rf
         #include "iwt_utils.c"        
   #undef TYPE_REAL     
      #define TYPE_DOUBLE
         #define TYPE_NAME 2cd
         #include "iwt_utils.c"
      #undef TYPE_DOUBLE
         #define TYPE_NAME 2cf
         #include "iwt_utils.c"           
        
        
                        
#endif
