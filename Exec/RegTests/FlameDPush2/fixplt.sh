for x in plt*
do
  if [ ${#x} = 9 ]
    then t="plt0"$(echo $x|sed -r 's/plt//')
    mv $x $t
    echo $t
  fi    
done

