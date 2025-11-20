      i <- 16
      f.num<-mp$f.numbers[i]
      feat = mp$features[,i]
      simplePlot(feat)
      feat.padded.ft.c = mp$features.padded.ft.c[,i]
      #
      refs = mp$refs
      refs.padded.ft = mp$refs.padded.ft
        r.num = 29
        ref = refs[,r.num, drop = F]
        ref.ft = refs.padded.ft[,r.num, drop = F]
        message(r.num)

                                  pad.size <- length(feat)-1
                                  feat.ft.c <- feat.padded.ft.c

