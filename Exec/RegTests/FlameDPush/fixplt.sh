for x in plt*
do
  if [ ${#x} = 8 ]
    then t="plt0"$(echo $x|sed -r 's/plt//')
    mv $x $t
    echo $t
  fi    
done

