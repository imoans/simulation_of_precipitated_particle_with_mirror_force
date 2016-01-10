#!/bin/bash

###
# 与えた平均と標準偏差を持つ正規分布からピッチ角をランダムに1000個生成し、
# そのピッチ角での衝突点を求める。
# nmf.d : ミラー力を考慮していない場合
# mf.d : ミラー力を考慮している場合
###

mean=$1
sd=$2

if [[ ! -n $mean ]]; then
    echo "arg1 is required: which is used for mean value of normal distribution of PA" 
    exit 1
fi

if [[ ! -n $sd ]]; then
    echo "arg2 is required: which is used for sd value of normal distribution of PA" 
    exit 1
fi

dir=data/collidedHeights/"$mean"_"$sd"

mkdir -p $dir

exe="./bin/collided_heights.out"

$exe --ignoreMF --useNorm > $dir/nmf.d &
pid1=$!
$exe            --useNorm > $dir/mf.d  &
pid2=$!

echo kill $pid1 $pid2
echo "watch wc $dir"
