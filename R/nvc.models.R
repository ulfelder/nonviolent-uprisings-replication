# CHENOWETH PROJECT: MODEL OBJECTS
# Jay Ulfelder
# 2013-12-21

gripes.f <- formula(nvc.start.1 ~ log(wdi.pop) +                        # Population size (log)
                    log(xxxcimr) +                                      # Infant mortality rate rel to global median (log)
                    wdi.gdpchg.s +                                      # GDP growth rate (sq rt)
                    sqrt(wdi.cpi) +                                     # Consumer price index (sq rt)
                    log1p(bnn.yroff) +                                  # Leader's years in office (log)
                    elceleth.c +                                        # Ruling elites' ethnicity is politically salient
                    dispota4.c +                                        # Any state-led discrimination
                    cir.physint +                                       # CIRI physical integrity index
                    I(cir.physint^2) )

modnzn.f <- formula(nvc.start.1 ~ log(wdi.pop) +                        # Population size (log)
                    wdi.popurb.mi +                                     # Urban population (0-100%), imputed
                    I(wdi.manuf.mi + wdi.servs.mi )+                    # Manufacturing & services (% of GDP), imputed
                    wdi.sch2.mi +                                       # 2ry school enrollment rate, imputed
                    log1p(wdi.mobp100) +                                # Mobile phone subs per 100 ppl (logged)
                    ios.gattwto )                                       # GATT/WTO member

rscmob.f <- formula(nvc.start.1 ~ log(wdi.pop) +                        # Population size (log)
                    wdi.popurb.mi +                                     # Urban population (0-100%), imputed
                    ythbul4 +                                           # Youth bulge (15-24 yos as % of tot pop)
                    log1p(bnk.unrest) +                                 # Sum of Banks riots & dems (log)
                    log1p(bnk.strikes) +                                # Sum of Banks strikes (log)
                    log1p(nvc.dosregt) +                                # Onsets of nonviolent campaigns in same region (log)
                    nvc.ongoing +                                       # Any ongoing nonviolent campaign
                    civilwar )                                          # Any ongoing PITF civil war

polopp.f <- formula(nvc.start.1 ~ log(wdi.pop) +                        # Population size (logged)
                    log1p(age) +                                        # Country age (logged)
                    postcoldwar +                                       # Post-Cold War period (1991+)
                    ios.iccpr1 +                                        # Signatory to ICCPR 1st Optional Protocol
                    nld.any.1 +                                         # Election year
                    pitfdem +                                           # PITF democracy indicator
                    I(pitfdem * nld.any.1) +                            # Election year/democracy interaction
                    I(postcoldwar * nld.any.1) +                        # Election year/post-cold war interaction
                    as.factor(fiw.cl) +                                 # FH civil liberties index (flipped & centered)
                    log1p(pol.durable) +                                # Regime durability
                    log1p(cou.tries5) )                                 # Coup activity in past 5 years

base.f <- formula(nvc.start.1 ~ log(wdi.pop))                           # Population size (log)
