################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../auxiliary.c \
../dispersal.c \
../epidemics.c \
../inout.c \
../param.c \
../pathogens.c \
../populations.c \
../sampling.c \
../sumstat.c 

OBJS += \
./auxiliary.o \
./dispersal.o \
./epidemics.o \
./inout.o \
./param.o \
./pathogens.o \
./populations.o \
./sampling.o \
./sumstat.o 

C_DEPS += \
./auxiliary.d \
./dispersal.d \
./epidemics.d \
./inout.d \
./param.d \
./pathogens.d \
./populations.d \
./sampling.d \
./sumstat.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


