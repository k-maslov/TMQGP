import QuarkTM

def get_channels(quark, gluon, T, G, G1, L, screen, suppress, lmax=1):
    params = {'G' : G, 'L' : L, 'screen' : screen}
    params1 = {'G' : G1, 'L' : L, 'screen' : screen}

    params_QG = {'G' : G, 'L': suppress*L, 'screen': screen}
    params_QG1 = {'G' : G1, 'L' : suppress*L, 'screen' : screen}

    params_GG = {'G' : G, 'L': suppress**2 * L, 'screen': screen}
    params_GG1 = {'G' : G1, 'L' : suppress**2 * L, 'screen' : screen}
    
    params_rep = params.copy()
    params_rep['G'] = -params_rep['G']
    params_rep1 = params1.copy()
    params_rep1['G'] = -params_rep1['G']
    params_rep_QG = params_QG.copy()
    params_rep_QG['G'] = -params_QG['G']
    params_rep_QG1 = params_QG1.copy()
    params_rep_QG1['G'] = -params_QG1['G']
    params_rep_GG = params_GG.copy()
    params_rep_GG['G'] = -params_GG['G']
    params_rep_GG1 = params_GG1.copy()
    params_rep_GG1['G'] = -params_GG1['G']
    channels_QQ = QuarkTM.ChannelGroup()
    channels_QG = QuarkTM.ChannelGroup()

    pss = [params, params1]
    pss_rep = [params_rep, params_rep1]

    pss_QG = [params_QG, params_QG1]
    pss_rep_QG = [params_rep_QG, params_rep_QG1]

    pss_GG = [params_GG, params_GG1]
    pss_rep_GG = [params_rep_GG, params_rep_GG1]

    labels = ['QQ', 'QG', 'GQ', 'GG']


    channels_QQ.addChannel(
        QuarkTM.ChannelL('qq3', lmax, quark, quark, T, pss, ds=4, da=3, Fa=1/2, )
    )

    channels_QQ.addChannel(
        QuarkTM.ChannelL('qa1', lmax, quark, quark, T, pss, ds=4, da=1, Fa=1, )
    )

    channels_QQ.addChannel(
        QuarkTM.ChannelL('qq6', lmax, quark, quark, T, pss_rep, ds=4, da=6, Fa=1/4, )
    )

    channels_QQ.addChannel(
        QuarkTM.ChannelL('qa8', lmax, quark, quark, T, pss_rep, ds=4, da=8, Fa=1/8, )
    )

    channels_QG.addChannel(
        QuarkTM.ChannelL('qg3', lmax, quark, gluon, T, pss_QG, ds=4, da=3, Fa=9./8, )
    )

    channels_QG.addChannel(
        QuarkTM.ChannelL('qg6', lmax, quark, gluon, T, pss_QG, ds=4, da=6, Fa=3./8, )
    )

    channels_QG.addChannel(
        QuarkTM.ChannelL('qg15', lmax, quark, gluon, T, pss_rep_QG, ds=4, da=15, Fa=3./8, )
    )

    channels_GQ = QuarkTM.ChannelGroup()
    channels_GG = QuarkTM.ChannelGroup()
# 
    channels_GQ.addChannel(
        QuarkTM.ChannelL('gq3', lmax, gluon, quark, T, pss_QG, ds=4, da=3, Fa=9/8)
    )

    channels_GQ.addChannel(
        QuarkTM.ChannelL('gq6', lmax, gluon, quark, T, pss_QG, ds=4, da=6, Fa=3/8)
    )

    channels_GQ.addChannel(
        QuarkTM.ChannelL('gq15', lmax, gluon, quark, T, pss_QG, ds=4, da=15, Fa=3/8)
    )

    channels_GG.addChannel(
        QuarkTM.ChannelL('gg1', lmax, gluon, gluon, T, pss_GG, ds=4, da=1, Fa=9/4)
    )

    channels_GG.addChannel(
        QuarkTM.ChannelL('gg16', lmax, gluon, gluon, T, pss_GG, ds=4, da=16, Fa=9/8)
    )

    channels_GG.addChannel(
        QuarkTM.ChannelL('gg27', lmax, gluon, gluon, T, pss_rep_GG, ds=4, da=27, Fa=3/4)
    )
