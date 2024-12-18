CC = gcc
CFLAGS = -Ilib

# ソースファイルのリスト
SRCS = ASPS_DSS.c lib/Common.c lib/ParameterDefined.c lib/Inverter_function.c 

# 出力ファイル
TARGET = output

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET)