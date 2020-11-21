
/* nd_setup helps with indexing the n-dim data
 * indices are organized into sets, with each set being one row/column/page/etc
 *
 * Travis Smith
 * Summer, 2011
 *
 * 
 * IN:
 * size = array containing the size of the data
 * ndims = length of the size array
 * opdim = dimension the operation will be performed in
 * OUT:
 * start = array of starting index of each set
 * nstart = length of start array
 * incr = increment between elements in the set
 */
void nd_setup(mwSize *size, mwSize ndims, mwSize opdim, mwSize **start, mwSize *nstart, mwSize *incr)
{
  mwSize ii, jj, kk, off;
  mwSize num_sets, num_in_each_set, inter_set_space;
  mwSize *start_ptr;
  
  inter_set_space = size[0];
  for (ii=1; ii<=opdim; ii++) { inter_set_space *= size[ii]; }
  
  num_sets = 1;
  for (ii=opdim+1; ii<ndims; ii++) { num_sets *= size[ii]; }  

  num_in_each_set = inter_set_space / size[opdim];
  *incr = num_in_each_set;
  
  *nstart = num_sets * num_in_each_set;
  *start = (mwSize *)mxMalloc(*nstart*sizeof(mwSize));
  start_ptr = *start;
  
  for (ii=0, kk=0; ii<num_sets; ii++)
  {
    off = ii*inter_set_space;
    for (jj=0; jj<num_in_each_set; jj++)
    {
      start_ptr[kk++] = off + jj;
    }
  }

  #ifdef VERBOSE
  printf("Size = %d",size[0]);
  for (ii=1; ii<ndims; ii++) { printf("x%d",size[ii]); }
  printf("\nopdim = %d\n",opdim);
  printf("nstart = %d, incr = %d\n",*nstart,*incr);
  #endif
}