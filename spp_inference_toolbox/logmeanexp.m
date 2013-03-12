function log_out = logmeanexp(x)

[y, maxx] = lowexp(x);

log_out = log(mean(y(:))) + maxx;



