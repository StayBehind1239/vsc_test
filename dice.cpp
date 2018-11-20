int roll(int i){

    if( i % 2 == 0)
    {
      i = (i*2)+1;
    }
    else
    {
      i++;
    }

    return i;
}