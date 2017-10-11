# Trying FLash 32 / 64 ADOLC with different OS
library(FLash)

data(ple4)
ple4_sr <- fmle(as.FLSR(ple4, model="bevholt"), control=list(trace=0))

# Mean F
f_status_quo <- 0.4

ctrl_target <- data.frame(year = 2000:2008,
              quantity = "f",
              val = f_status_quo)
ctrl_f <- fwdControl(ctrl_target)

ple4_f_sq <- fwd(ple4, ctrl = ctrl_f, sr = ple4_sr)
plot(ple4_f_sq)
fbar(ple4_f_sq)

# Linux Orig - looks OK

# Constant catch 1
future_catch <- 1000
ctrl_catch <- fwdControl(
    data.frame(
        year=2000:2008,
        quantity = "catch",
        val=future_catch))

ple4_catch <- fwd(ple4, ctrl_catch, sr = ple4_sr, maxF=5)
plot(ple4_catch)
catch(ple4_catch)

# Linux Orig - Is OK
# Linux new - Is OK
# Windows 64 - Is OK

# Constant catch 2
future_catch <- 10000
ctrl_catch <- fwdControl(
    data.frame(
        year=2000:2008,
        quantity = "catch",
        val=future_catch))

ple4_catch <- fwd(ple4, ctrl_catch, sr = ple4_sr, maxF=5)
plot(ple4_catch)
catch(ple4_catch)
fbar(ple4_catch)

# Linux Orig - is OK
# Linux New - is OK
# Windows 64 - is OK

# Constant catch 3
future_catch <- c(catch(ple4)[,ac(2000:2008)])
ctrl_catch <- fwdControl(
    data.frame(
        year=2000:2008,
        quantity = "catch",
        val=future_catch))

ple4_catch <- fwd(ple4, ctrl_catch, sr = ple4_sr, maxF=5)
plot(ple4_catch)
catch(ple4_catch)[,ac(2000:2008)] - future_catch
fbar(ple4_catch)

# Linux Orig - is OK
# Linux New - is OK
# Windows 64 - is OK


