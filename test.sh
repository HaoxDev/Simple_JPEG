#QF = 5
echo "QF:5"
./encode.o ./Baboon.raw ./encode_tmp.o 5
./decode.o ./encode_tmp.o ./tryresult.raw 5
python ./psnr.py ./tryresult.raw ./Baboon.raw

#QF = 10
echo "QF:10"
./encode.o ./Baboon.raw ./encode_tmp.o 10
./decode.o ./encode_tmp.o ./tryresult.raw 10
python ./psnr.py ./tryresult.raw ./Baboon.raw

#QF = 20
echo "QF:20"
./encode.o ./Baboon.raw ./encode_tmp.o 20
./decode.o ./encode_tmp.o ./tryresult.raw 20
python ./psnr.py ./tryresult.raw ./Baboon.raw

#QF = 50
echo "QF:50"
./encode.o ./Baboon.raw ./encode_tmp.o 50
./decode.o ./encode_tmp.o ./tryresult.raw 50
python ./psnr.py ./tryresult.raw ./Baboon.raw

#QF = 80
echo "QF:80"
./encode.o ./Baboon.raw ./encode_tmp.o 80
./decode.o ./encode_tmp.o ./tryresult.raw 80
python ./psnr.py ./tryresult.raw ./Baboon.raw

#QF = 90
echo "QF:90"
./encode.o ./Baboon.raw ./encode_tmp.o 90
./decode.o ./encode_tmp.o ./tryresult.raw 90
python ./psnr.py ./tryresult.raw ./Baboon.raw
