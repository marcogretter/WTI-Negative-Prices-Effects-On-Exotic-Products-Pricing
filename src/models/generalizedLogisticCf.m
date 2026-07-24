function phi = generalizedLogisticCf(u, alpha, betaPar)

logPhi = generalizedLogisticLogCf(u, alpha, betaPar);

phi = exp(logPhi);
end