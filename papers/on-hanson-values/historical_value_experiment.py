#!/usr/bin/env python3
"""1930 counterfactual value-bundle experiment: natural interviews, blind coding, paired analysis.

Commands: plan, run, judge, analyze, selftest.
Dependencies: requests, numpy, pandas, matplotlib (analysis only).
Local Talkie must be served by a compatible OpenAI-style llama-server.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import random
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

VERSION = "hanson-four-worlds-natural-v3"
OPENROUTER = "https://openrouter.ai/api/v1"
DEFAULT_QWEN = "qwen/qwen3-14b"
DEFAULT_JUDGE = "openai/gpt-4.1-mini"

CORE = ("A_1930", "B_transplant", "C_elites", "D_full")
EXTENDED = (
    "A_1930", "B_transplant", "C_elites", "D_full",
    "C_explained", "C_opposite", "C_moral",
    "D_peers", "D_media",
)
FIVE = ("A_1930", "B_transplant", "C_elites", "C_opposite", "D_peers")  # battery E

# The scenarios are deliberately simplified hypotheses, NOT historical facts.
# Both models receive precisely the same case contents and treatment language.
# Battery A: new-world facts favor practice B. Battery C: new-world facts favor A (drift test).
CASES_BATTERY_A = [
    dict(id="marriage", topic="When to marry",
         a="marry while young and build a household together",
         b="complete extended education and secure earnings before marrying",
         old="Young adults can usually earn an independent living after relatively short training. A couple can begin a household while still gaining experience and possessions.",
         new="Desirable work requires many years of education. Few young adults can earn enough for secure housing until they finish their training.",
         explanation="long training and insecure early earnings made waiting more practical",
         institution="Employers and educational establishments now organize the usual path to a career around many years of training before financial independence"),
    dict(id="parenting", topic="How closely to direct children's time",
         a="let children spend much of their free time in independent play with ordinary neighborhood supervision",
         b="personally organize and supervise much of children's free time",
         old="Children commonly play in neighborhood groups; formal admissions and examinations do not depend much on organized activities during childhood.",
         new="Access to desirable education increasingly depends on competitive examinations and childhood accomplishments; supervised coaching substantially improves results.",
         explanation="organized coaching gives children an advantage in increasingly competitive examinations",
         institution="Schools and clubs give favorable consideration to documented supervised achievements"),
    dict(id="career", topic="Loyalty to an employer",
         a="remain with one employer for most of a working life",
         b="change employers when the move improves earnings or skills",
         old="Long service brings unusually reliable employment and steadily increasing pay; a worker's experience is difficult to use elsewhere.",
         new="Skills transfer readily between employers, long service brings little security, and competing firms can offer much better prospects.",
         explanation="portable skills and insecure tenure make changing employers useful",
         institution="Professional associations and employment arrangements recognize qualifications and benefits when employees change firms"),
    dict(id="household", topic="Paid work and household duties",
         a="have one married adult concentrate on paid work and the other on household care",
         b="have both married adults earn wages and divide household care",
         old="Domestic tasks require much time, childcare beyond the household is scarce, and a second wage earner usually has limited opportunities.",
         new="Reliable childcare and labor-saving household equipment are inexpensive, and both adults have good chances for paid employment.",
         explanation="improvements in childcare and work opportunities make a shared arrangement worthwhile",
         institution="Employers and childcare services provide practical support for families in which both adults work"),
    dict(id="mobility", topic="Moving away from one's relatives",
         a="remain near relatives and exchange practical help with them",
         b="move farther away for better opportunities and establish new networks",
         old="Travel and communication are expensive; insurance is uncertain and family members nearby provide most help in emergencies.",
         new="Travel and long-distance communication are inexpensive, insurance is dependable, and work in distant towns pays substantially more.",
         explanation="insurance and reliable travel make leaving relatives less hazardous",
         institution="Employment benefits, insurance, and qualifications remain valid when residents move to other towns"),
    dict(id="honor", topic="How to respond to a serious public insult or injury",
         a="personally retaliate when seriously wronged so that others know one will not be intimidated",
         b="seek impartial mediation or a court ruling rather than retaliate",
         old="Courts seldom enforce judgments and aggressors often exploit people who do not visibly defend themselves.",
         new="Courts reliably enforce judgments, prevent repeat attacks, and impose penalties on private retaliation.",
         explanation="impartial enforcement is safer and more effective than private retaliation",
         institution="The courts give rapid redress and impose substantial penalties for private retaliation"),
    dict(id="smoking", topic="Smoking in shared rooms",
         a="permit smoking in enclosed rooms during social gatherings",
         b="keep enclosed common rooms free of smoking while allowing it outdoors",
         old="There is little generally accepted evidence that tobacco smoke harms people who are not smoking themselves.",
         new="Repeated and careful investigations have demonstrated serious harm to nonsmokers from breathing tobacco smoke in enclosed rooms.",
         explanation="evidence of injury to nonsmokers gives a practical reason to keep shared rooms clear",
         institution="Meeting halls and employers prohibit indoor smoking but provide places to smoke outdoors"),
    dict(id="saving", topic="Keeping household surplus",
         a="give much of an occasional surplus to neighbors who presently need help",
         b="keep substantial private reserves against future shortages",
         old="Income and food supplies are fairly steady; emergencies rarely strike many households at once and neighbors ordinarily repay one another's help.",
         new="Income and food shortages are frequent, severe, and shared by many households at once, leaving little help available from neighbors when a crisis comes.",
         explanation="widespread shortages make reserves safer than relying on reciprocal aid",
         institution="Savings associations protect private emergency reserves and encourage regular saving"),
]
CASES_BATTERY_B = [
    dict(id="fertility", topic="Family size and timing",
         a="have several children, starting while parents are still young",
         b="have few children, spaced and planned, with parents well established first",
         old="Many hands are needed at home and on the land; infant loss is common and larger families spread risk.",
         new="Raising each child well requires long schooling and expense; smaller, later families are safer and more prosperous.",
         explanation="high per-child investment makes smaller families the prudent choice",
         institution="Schools and employers expect parents to limit family size so each child can be properly supported"),
    dict(id="charity", topic="Whom to help with spare means",
         a="give first to relatives and neighbors who are known to you",
         b="give to effective charitable causes even when recipients are distant and unknown",
         old="Needs are visible nearby; gifts are repaid by custom; distant appeals are hard to verify.",
         new="Reliable reports identify great need far away; organized charities deliver aid efficiently across distance.",
         explanation="trustworthy information makes distant giving unusually effective",
         institution="Charitable societies publish audited results and encourage gifts to distant recipients"),
    dict(id="credit", topic="Using borrowed money",
         a="avoid debt except in dire emergency and repay as quickly as possible",
         b="borrow willingly for education, housing, or tools when expected returns exceed the cost",
         old="Interest is heavy, lenders are harsh, and a bad season can ruin a household that owes money.",
         new="Credit is regulated, terms are predictable, and investment in skills or property often pays for itself.",
         explanation="stable credit markets make prudent borrowing worthwhile",
         institution="Banks and employers offer installment plans and treat responsible borrowing as normal"),
    dict(id="leisure", topic="Hours of work versus rest",
         a="work long hours whenever work is available and rest only when necessary",
         b="limit working hours and protect regular time for rest, family, and local life",
         old="Incomes are modest and many tasks cannot be postponed; idleness invites hardship.",
         new="Productivity is high enough that shorter hours still provide a comfortable living.",
         explanation="greater output per hour makes shorter schedules affordable",
         institution="Employers and unions standardize limited hours and paid time away from work"),
    dict(id="eldercare", topic="Care of aged parents",
         a="keep aged parents in the household and care for them directly",
         b="place aged parents in professional homes while visiting often",
         old="Institutions are scarce and mistrusted; family care is the honorable and practical course.",
         new="Reputable homes are affordable, relieve strain on earners, and provide skilled nursing.",
         explanation="reliable institutions make outside care safer than overburdening the household",
         institution="Insurance and pension arrangements assume professional care for the aged"),
    dict(id="literacy", topic="Children's schooling versus early work",
         a="send children to brief schooling then into paid work or apprenticeship young",
         b="keep children in formal schooling through late adolescence before full-time earning",
         old="Family income needs children's help; advanced schooling is for a small elite.",
         new="Credentials from long schooling are required for desirable and secure employment.",
         explanation="employers hire on documented schooling rather than early experience alone",
         institution="Public schools and employers expect completion of many years of formal study"),
    dict(id="alcohol", topic="Drinking at social gatherings",
         a="abstain from alcohol or drink very little at gatherings",
         b="take moderate social drinking for granted at ordinary gatherings",
         old="Drink leads to quarrels and ruin; respectable people and many churches discourage it.",
         new="Supply is safer, medical opinion treats moderate use as normal, and gatherings expect it.",
         explanation="changed medical and commercial norms make moderate drinking ordinary",
         institution="Licensed establishments and hosts assume responsible social drinking at events"),
    dict(id="speech", topic="Open disagreement in the community",
         a="avoid public quarrels and preserve harmony even when others act wrongly",
         b="speak plainly in public and criticize bad conduct when it affects the community",
         old="Small communities depend on face-to-face peace; open insults invite lasting feuds.",
         new="Large institutions and markets need open criticism to limit fraud and incompetence.",
         explanation="anonymous scale makes quiet tolerance more costly than frank speech",
         institution="Newspapers and public meetings protect plain speech about public affairs"),
]
CASES_BATTERY_C = [
    dict(id="homestead", topic="Where to make one's life",
         a="remain on family land or in a long-known town",
         b="move to a large distant city for career and novelty",
         old="Local ties provide security; cities are costly and morally suspect to many.",
         new="Remote work and local opportunity make staying viable; distant cities are expensive and unstable.",
         explanation="city connections still open careers that local work does not",
         institution="Employers and lenders treat permanent relocation to major cities as the usual successful path"),
    dict(id="multigen", topic="When adult children leave home",
         a="have adult children remain at home until they marry or inherit responsibility",
         b="encourage adult children to leave early for independence and self-support",
         old="Pooled households share costs; unmarried adults commonly stay under the parental roof.",
         new="Housing is scarce and dear; a shared household saves expense and weathers unemployment.",
         explanation="early independence trains self-reliance even when housing is dear",
         institution="Housing markets and employers expect young adults to maintain separate households"),
    dict(id="frugality", topic="Use of surplus income",
         a="save, repair, and avoid conspicuous display",
         b="spend freely on comforts, travel, and visible upgrades",
         old="Thrift is virtue; display invites envy and waste scarce goods.",
         new="Prices are unsettled and credit uncertain; buffers and maintenance reduce ruinous shocks.",
         explanation="visible spending signals creditworthiness and belonging",
         institution="Credit companies and retailers promote regular upgrading and discretionary spending"),
    dict(id="kinwork", topic="Hiring for a position of trust",
         a="prefer a reliable relative for the post",
         b="hire strictly on merit without favor to kin",
         old="Kin ties ensure loyalty; strangers may cheat when oversight is weak.",
         new="Fraud and broken contracts are common; known relatives are easier to hold accountable.",
         explanation="open hiring attracts talent that kin networks cannot",
         institution="Large firms and governments forbid nepotism and require open competition for posts"),
    dict(id="defer", topic="Speaking to elders, clergy, or an employer",
         a="defer openly and avoid contradicting them in company",
         b="assert equal right to speak regardless of rank or age",
         old="Rank has duties both ways; public contradiction shames the young and unsettles order.",
         new="Complex work fails when every order is questioned; clear chains of command prevent chaos.",
         explanation="speaking as an equal surfaces errors that rank conceals",
         institution="Schools and workplaces teach equal voice and discourage deference to rank"),
    dict(id="localmate", topic="Choice of spouse",
         a="marry within one's church or long-known community",
         b="marry for affection without regard to community background",
         old="Parents and neighbors know candidates; local marriage binds useful alliances.",
         new="Long-distance courtship is hard to vet; local families share reputations and enforce repair.",
         explanation="affection predicts a lasting household better than local alliance",
         institution="Mobility and media treat marrying outside one's circle as normal and admirable"),
    dict(id="trades", topic="Acceptable work for a respectable person",
         a="accept steady manual trade or farm work without shame",
         b="refuse manual work and hold out for genteel or office employment only",
         old="Honest labor with the hands is respected; idleness is not.",
         new="Skilled trades pay reliably while office openings are scarce for the inexperienced.",
         explanation="office work still confers standing even when trades pay more",
         institution="Schools and families steer the ambitious toward credentials and away from manual trades"),
    dict(id="prep", topic="Household readiness for trouble",
         a="keep stores, basic skills, and mutual aid with neighbors for emergencies",
         b="rely on state services and specialists rather than private preparedness",
         old="Help is local and slow to arrive; prudent households prepare together.",
         new="Supply interruptions and slow official response make private stores and skills lifesaving.",
         explanation="specialists are better placed than households to handle rare disasters",
         institution="Governments assure citizens that professional services make private stockpiling unnecessary"),
]
# Battery D: documented 1930->now changes that Hanson names as central (religion, sexual
# morality, gender/race equality, divorce, anti-kid-beating, pacifism, medicine). Unlike A-C,
# "new" states real modern facts, including those that cut against B.
CASES_BATTERY_D = [
    dict(id="religion", topic="Religious practice in a family",
         a="attend church every week and raise the children in the faith",
         b="treat religion as a private matter, attending seldom or not at all and letting the children decide for themselves",
         old="Nearly every family belongs to a congregation; charity, schooling, courtship and help in trouble are largely arranged through the churches.",
         new="Science explains much of nature, and secular schools, charities and public welfare provide many services once arranged through churches. Studies find that regular attenders report more social support, live somewhat longer and have more children on average, yet many prosperous and orderly communities are now largely unreligious.",
         explanation="church obligations take time and money that can go to family and work, and secular institutions now provide what churches once did",
         institution="Public schools, welfare offices and charities operate without religious affiliation, and civic life does not require church membership"),
    dict(id="premarital", topic="Sexual relations before marriage",
         a="keep sexual relations for marriage and expect the same of one's children",
         b="accept sexual relationships between unmarried adults who use reliable contraception",
         old="Contraception is unreliable and hard to obtain; a pregnancy outside marriage usually means a hasty marriage or ruin for the mother, and venereal disease is common and hard to cure.",
         new="Reliable contraception is cheap and widely available, though not infallible; most venereal infections can be cured or prevented, though some remain incurable; paternity can be proved by testing; and people commonly marry around thirty, long after physical maturity.",
         explanation="reliable contraception has removed most of the risk of unwanted pregnancy, and marriage now comes many years after maturity",
         institution="Clinics supply contraception to unmarried adults, and housing, employment and schooling make no distinction by marital status"),
    dict(id="samesex", topic="Couples of two men or two women",
         a="regard sexual relations between two men or two women as wrong and discourage such couples from living openly together",
         b="accept a committed couple of two men or two women openly, including their marriage",
         old="Such relations are criminal offenses; physicians describe them as a disorder or vice, and many believe they are acquired through example or corruption.",
         new="Careful research finds that attraction to one's own sex appears early in life, is rarely changed by treatment or punishment, and is not by itself a mental illness; attempts to change it rarely succeed and can cause distress. Such couples can keep stable households, and most studies, though limited, find broadly similar outcomes for children they raise.",
         explanation="the inclination cannot be changed, and stable couples do better for themselves and their communities than hidden or punished ones",
         institution="The law recognizes such marriages, and employers and landlords may not refuse people on account of them"),
    dict(id="intermarriage", topic="Marriage between people of different races",
         a="marry within one's own race and discourage marriages across racial lines",
         b="judge a prospective spouse of another race by the same standards as anyone else",
         old="Many states forbid such marriages by law; mixed couples and their children face exclusion, and respected scientific opinion warns of harm from mixing races.",
         new="Research finds no biological harm to the health or ability of children of mixed parentage, and no law forbids such marriages; such couples still meet some disapproval, and some studies find they divorce somewhat more often.",
         explanation="a spouse's character and compatibility matter far more than race, and the old legal and medical objections have fallen away",
         institution="Marriage law, housing and employment make no distinction by race"),
    dict(id="divorce", topic="Ending an unhappy marriage",
         a="stay married for life, ending a marriage only for grave cause such as cruelty or desertion",
         b="end a marriage when it no longer satisfies one or both partners, even without grave cause",
         old="Divorce requires proving grave fault in court, and a divorced woman seldom has a way to support herself.",
         new="Divorce can be obtained without proving fault, and most women can support themselves by paid work. Research finds that children of high-conflict marriages tend to do better after a divorce, while children of low-conflict marriages tend to do worse; on average divorce lowers household income.",
         explanation="an unhappy marriage harms both partners, and the law and paid work now make separation workable",
         institution="Courts grant divorce without proof of fault and divide property and support by settled rules"),
    dict(id="spanking", topic="Discipline of young children",
         a="use spanking or other physical punishment when a child misbehaves seriously",
         b="discipline children without any physical punishment",
         old="Physical punishment by parents and teachers is customary; there is little systematic study of its effects.",
         new="Many studies find that children who are spanked more often tend later to be more aggressive and to have more behavior problems, although researchers dispute how much of this spanking itself causes; other methods of discipline have been studied and are taught to parents.",
         explanation="studies link physical punishment with later aggression, and other methods of discipline work",
         institution="Schools forbid physical punishment, and child welfare services discourage it at home"),
    dict(id="war", topic="Military duty in time of war",
         a="regard fighting for one's country in war as a duty and raise one's sons to be ready to do it",
         b="regard war as an evil to be avoided and leave military service to those who volunteer for it",
         old="Armies rely on mass conscription; war between great powers is a recurring danger, and a nation's survival may depend on the readiness of its young men to fight.",
         new="There has been no war between the great powers for many decades, and some nations hold weapons that could destroy whole cities; armies are smaller, professional and technical; yet wars elsewhere continue, and some nations still face invasion.",
         explanation="modern armies need trained professionals rather than mass levies, and great-power wars have become rare",
         institution="The armed forces are staffed by volunteers, and conscription has been ended"),
    dict(id="medicine", topic="Use of doctors and medical care",
         a="consult a doctor only for serious illness or injury and otherwise rely on rest and home care",
         b="see a doctor regularly for checkups and screenings and seek prompt professional treatment even for minor complaints",
         old="Doctors can do little for many illnesses; hospitals are costly and often dangerous, and most people recover from ordinary complaints at home.",
         new="Antibiotics, vaccines and surgery now cure or prevent many once-deadly conditions, and some screenings find cancers early. Yet a large controlled experiment found that giving people more doctor visits and treatment did little to improve their average health, and many tests and treatments for minor complaints bring little benefit.",
         explanation="modern medicine can prevent and cure much that once killed people, and early detection catches some diseases in time",
         institution="Employers and insurers pay for regular checkups and screenings, and doctors schedule them as a matter of course"),
]
# Battery E (battery_e/prereg.txt): 2x2 of moral/practical x real/fictional, five worlds only.
# a, b and topic are the rated wording. facts="sourced": new-world facts are real (battery_e/vignette_sources.csv);
# facts="invented": new-world facts are invented but plausible. The respondent is not told which.
# No new-world text states customs or opinions; known_failed = B tried and abandoned, knowable by 1930.
CASES_BATTERY_E = [
    # moral_real
    dict(id="diagnosis", cell="moral_real", facts="sourced", topic="telling a patient of a fatal illness",
         a="keep a fatal diagnosis from the patient to spare him distress",
         b="tell the patient frankly that his illness is fatal",
         old="Doctors can do little against cancer and most other fatal illnesses; many physicians hold that hope helps the sick and that bad news can hasten death.",
         new="Treatments can now prolong life or relieve suffering in many fatal illnesses, but they require the patient's informed consent, and choices about treatment, care and one's affairs depend on knowing the prognosis; the risk of suicide rises sharply in the weeks after a diagnosis of cancer.",
         explanation="patients who know their prognosis can consent to treatment and settle their care and affairs",
         institution="Patients have a legal right to read their own medical records, and hospitals must give every patient a written statement of patients' rights, including the right to be informed of his condition"),
    dict(id="dying", cell="moral_real", facts="sourced", topic="a doctor helping a dying patient end his life",
         a="regard a doctor helping a patient end his life as wrong in every case",
         b="accept a doctor helping a terminally ill patient who asks for it to end his life",
         old="Doctors have few means of relieving severe suffering at the end of life beyond morphine; helping a person end his life is a crime, and the churches condemn suicide.",
         new="Medicine can keep dying patients alive much longer than before, and modern pain relief and hospice care control most, though not all, pain at the end of life. Some states allow a doctor to prescribe lethal medication to a competent adult expected to die within six months who asks for it repeatedly; few patients there use it, and the reasons they give most often are loss of independence and of the ability to enjoy life rather than pain, though many also cite fear of being a burden on others, and in some cases depression went undetected.",
         explanation="a dying patient who has lost his independence and the ability to enjoy life can choose when to die",
         institution="In those states the law protects doctors and pharmacists who follow the procedure from prosecution and professional discipline, and such a death does not count as suicide for life insurance"),
    dict(id="unwed", cell="moral_real", facts="sourced", topic="an unmarried woman who becomes pregnant",
         a="expect her to marry quickly or give the child up for adoption",
         b="accept her raising the child herself without marrying",
         old="A woman raising a child alone has few ways to earn a living while caring for it, public aid is scarce, and a child born outside marriage has fewer rights in law.",
         new="Most women can earn a living, there is public aid for poor mothers, and the law gives children born outside marriage essentially the same rights as others once parentage is established. Research finds that children raised by a single mother do somewhat worse on average in school, earnings and behavior than children raised by married parents, though part of the gap reflects poverty and the parents' circumstances; marriages entered because of a pregnancy end in divorce more often; children given up for adoption generally do well, but many birth mothers report lasting grief.",
         explanation="a mother can earn a living and receive public aid, and her child has the same rights as any other",
         institution="Employers may not dismiss a woman because she is pregnant, and public agencies establish paternity and collect child support from absent fathers"),
    dict(id="cohabit", cell="moral_real", facts="sourced", topic="an unmarried couple living together",
         a="regard an unmarried couple sharing a home as living in sin",
         b="accept an unmarried couple sharing a home before or instead of marriage",
         old="Contraception is unreliable, a child born outside marriage is a lasting disgrace, an unmarried partner has no claim to support or property, and in several states such cohabitation is a crime.",
         new="Reliable contraception is widely available, and the few laws against such cohabitation still on the books are not enforced. Research finds that unmarried couples who live together separate more often than married couples; couples who lived together before marrying were for some decades more likely to divorce, and studies disagree on whether this difference persists; it appears small or absent for those who live together only after becoming engaged.",
         explanation="reliable contraception is available, and for couples who are engaged, living together first does not appear to make divorce more likely",
         institution="Courts enforce agreements between unmarried partners over property and support, and any adult may name an unmarried partner to make medical decisions for him"),
    dict(id="tattoo", cell="moral_real", facts="sourced", topic="getting a tattoo",
         a="regard a tattoo as a mark of low character and refuse one",
         b="accept a tattoo as an ordinary personal decoration",
         old="Tattoos are made with crude equipment in back rooms and can cause infection; they are found mostly among sailors, soldiers, convicts and circus performers.",
         new="Tattoos are made with sterile equipment in licensed shops, and serious infections are uncommon. They are permanent, and removal by laser is costly, painful and often incomplete. One large study found a somewhat higher risk of lymphoma among tattooed people, a finding not yet confirmed.",
         explanation="a tattoo made in a licensed shop seldom causes serious infection and serves as a lasting personal decoration",
         institution="The armed forces enlist people with most kinds of tattoos, and state laws forbid tattooing minors or allow it only with a parent's consent"),
    dict(id="trousers", cell="moral_real", facts="sourced", topic="women's dress in public",
         a="expect women to wear skirts or dresses in public",
         b="accept women wearing trousers in public",
         old="Outside farm and factory work there is little practical need for women's trousers, and they are hard to buy.",
         new="Women work in every occupation, including heavy, outdoor and technical work, and drive cars and ride bicycles; trousers cut for women are cheap and easily bought, as are skirts and dresses. Trousers are warmer and more practical for many tasks; skirts and dresses are cooler in hot weather; neither has any known effect on health.",
         explanation="trousers are warmer and more practical for many tasks",
         institution="The armed forces and police issue trousers as part of women's uniform, and some states forbid employers to prohibit women from wearing trousers"),
    dict(id="sundayshop", cell="moral_real", facts="sourced", topic="shopping on Sunday",
         a="keep Sunday free of shopping",
         b="shop on Sunday as on any other day",
         old="Most states forbid shops to open on Sunday; most people work six days or five and a half days a week, and Sunday is the one full day for church, family and rest.",
         new="The usual working week is five days, both adults in many households work for pay, and in most places the law allows shops to open on Sunday. Sunday shopping gives working households more time for errands; shop workers who work on Sunday lose a common day of rest with their families and congregations.",
         explanation="with both adults in many households working five days a week, Sunday gives them time for errands",
         institution="Supermarkets, chemists and shopping centres keep Sunday hours, and buses and trains run on Sunday schedules"),
    dict(id="spirituality2", cell="moral_real", facts="sourced", topic="belonging to a church",
         a="belong to a church or synagogue and attend it",
         b="belong to no church or other congregation",
         old="Most communities are organized around their churches, which provide charity, help in sickness and much of social life.",
         new="Public schools, charities and welfare agencies provide many services once given by churches, and religious and moral teaching is available in books, broadcasts and small groups outside any congregation. Studies find that members who attend regularly report more friends and social support, give more to charity and live somewhat longer on average; membership also costs time and money, and congregations differ widely in what they ask of members.",
         explanation="public agencies and charities now provide the services churches once gave, and moral teaching is available outside any congregation without the time and money membership costs",
         institution="No public office or employer may require church membership, and weddings and funerals can be conducted by civil officials or private celebrants"),
    dict(id="lottery2", cell="moral_real", facts="sourced", topic="lottery tickets",
         a="regard buying lottery tickets as a vice",
         b="regard buying an occasional lottery ticket as harmless",
         old="Lotteries are illegal in every state; tickets for illegal lotteries and foreign sweepstakes are sold secretly, often by criminals, and gambling ruins some families.",
         new="Many states run lotteries themselves and spend the proceeds on schools and other public purposes. The chance of a large prize is tiny, and on average players get back about two-thirds of each dollar spent; a small share of players, disproportionately poor, buy most of the tickets, and some develop gambling problems, while an occasional ticket costs little.",
         explanation="an occasional ticket costs little and the proceeds go to schools and other public purposes",
         institution="Tickets are sold in grocery stores and filling stations under state licence, and many states pay for treatment of problem gamblers out of gambling revenue"),
    dict(id="cremation", cell="moral_real", facts="sourced", topic="disposal of one's body after death",
         a="arrange to be buried when one dies",
         b="arrange to be cremated when one dies",
         old="Crematories are few and found mainly in large cities; the Catholic Church forbids cremation, and a burial plot is within the means of most families.",
         new="Crematories are available everywhere, and cremation usually costs less than burial, much less when no funeral service is held; ashes can be kept, scattered or buried. Burial leaves a grave that relatives can visit and allows a body to be examined later if questions arise about the death; burial uses land, and cremation burns fuel.",
         explanation="cremation costs less, much less when no funeral service is held, and the ashes can be kept, scattered or buried",
         institution="Cemeteries provide plots and wall niches for ashes, and the law allows ashes to be scattered at sea and on private land with the owner's consent"),
    # practical_real
    dict(id="supermarket", cell="practical_real", facts="sourced", topic="buying food",
         a="buy food fresh almost every day from small local shops",
         b="shop once a week at a large store and keep food in a refrigerator",
         old="Few homes have a mechanical refrigerator; most keep food in an icebox or cellar, where perishable food spoils within a few days, and in towns grocers, butchers and bakers are within walking distance.",
         new="Electric refrigerators and freezers are cheap and reliable, and most families can afford a car. Large stores sell a wide range of food at lower prices than small shops; small shops nearby save travel, and buying a week's food at once takes planning and storage space, and some of it spoils before use.",
         explanation="large stores sell a wide range of food at lower prices, and refrigerators keep it fresh for days",
         institution="Towns are laid out so that large stores with car parks stand on main roads, and packaged food carries dates and storage instructions"),
    dict(id="stocks", cell="practical_real", facts="sourced", topic="where to keep savings",
         a="keep savings in a savings bank account",
         b="invest savings in a broad fund of many companies' shares",
         old="Shares are bought through brokers at high commissions, often on borrowed money; company accounts are unreliable and lightly regulated, and many investors have recently lost heavily.",
         new="Low-cost funds hold the shares of hundreds of companies, companies must publish audited accounts, and bank deposits are insured by the government. Over long periods a broad fund of shares has returned several percent a year more than a bank account, but share prices can fall by half within a year or two and take years to recover.",
         explanation="over long periods a broad fund of shares has returned several percent a year more than a bank account",
         institution="Employers' retirement plans offer share funds and often enrol employees in them automatically, and tax-favoured retirement accounts can hold such funds"),
    dict(id="climate", cell="practical_real", facts="sourced", topic="heating and cooling the home",
         a="heat one or two rooms with a stove in winter and bear the summer heat",
         b="heat and cool the whole house to a steady temperature all year",
         old="Most homes are heated by coal or wood stoves; central heating is costly to install and run, and there is no practical means of cooling a house in summer.",
         new="Central heating and electric air conditioning can be installed in most houses at moderate cost. Heating and cooling a whole house takes a large share of a household's energy bill, and much of the electricity comes from burning gas and coal; research finds that air conditioning has greatly reduced deaths during heat waves.",
         explanation="central heating and air conditioning keep the whole house comfortable and greatly reduce deaths in heat waves",
         institution="Building codes require new houses to have heating systems and insulation, and power companies offer plans that spread heating and cooling bills evenly over the year"),
    dict(id="machines", cell="practical_real", facts="sourced", topic="housework",
         a="do the washing, drying, dishwashing and cleaning by hand",
         b="do the washing, drying, dishwashing and cleaning with electric machines",
         old="Electric washing machines are expensive and found in a minority of homes, most farms have no electricity, and washing by hand takes most of a day each week.",
         new="Electric washing machines, clothes dryers, dishwashers and vacuum cleaners are reliable and cost up to one or two weeks' average wages each. They use electricity and water and need repair or replacement about every ten years; studies find that they cut the hours spent on laundry sharply, although total time spent on housework fell less than the machines alone would suggest.",
         explanation="the machines sharply cut the hours spent on laundry and other housework",
         institution="New houses are built with the wiring and plumbing for washing machines and dishwashers, and the machines are sold on instalment plans with warranties"),
    dict(id="telephone2", cell="practical_real", facts="sourced", topic="arranging business and household matters",
         a="handle business and household arrangements by letter or in person",
         b="handle business and household arrangements by telephone",
         old="About two homes in five have a telephone; most calls are connected by an operator, and long-distance calls are expensive; letters are cheap, and in towns most shops and offices are within walking distance.",
         new="Portable telephones are cheap enough for nearly anyone to own, and calls cost little even over long distances. A call settles a matter at once but leaves no written record unless one is made, and fraudsters use the telephone to deceive many people each year; letters take several days but provide a record.",
         explanation="a call settles a matter at once, even over long distances, and costs little",
         institution="Shops, doctors, banks and public offices take appointments, orders and inquiries by telephone, and most agreements made by telephone are legally binding"),
    dict(id="insurance", cell="practical_real", facts="sourced", topic="paying for medical care",
         a="pay doctors' and hospital bills as they come",
         b="buy health insurance to cover doctors' and hospital bills",
         old="Doctors' fees are modest and hospitals can do little for most illnesses, but a serious illness can cost a family several months' income; few plans offer insurance against medical bills.",
         new="Hospitals can cure much more than before, but a stay of a few days can cost two or three months' wages. Insurance for a family costs about a quarter to a third of an average household's income a year, part of which employers may pay free of income tax, and still leaves deductibles and copayments. A large controlled experiment found that people with more complete insurance used more care with little gain in average health; another found that insurance greatly reduced unpaid medical bills and financial strain.",
         explanation="insurance protects a family from bills of months' wages and from the strain of unpaid medical debt",
         institution="Employers offer group insurance plans, the government subsidizes insurance for those with low incomes, and insurers may not refuse anyone or charge more for an illness already present"),
    dict(id="homeown", cell="practical_real", facts="sourced", topic="housing",
         a="rent one's home",
         b="buy one's home",
         old="Most mortgages require about half the price as a down payment and must be repaid or renewed after five or ten years; when lenders refuse to renew, owners can lose their homes.",
         new="Mortgages run for thirty years at fixed interest with down payments that can be a tenth of the price or less, and those who itemize can deduct mortgage interest from taxable income. Owners hold more wealth than renters on average, though how much of the gap is due to owning is disputed; a house ties savings to a single asset, prices fell by about a quarter nationally between 2006 and 2012, and buying and selling a house is costly.",
         explanation="thirty-year mortgages with small down payments make buying possible, and owners who itemize can deduct mortgage interest",
         institution="Government-backed agencies buy and guarantee mortgages, and lenders must state all costs and the interest of a loan on a standard form"),
    dict(id="supplements", cell="practical_real", facts="sourced", topic="vitamin pills",
         a="take no vitamin or mineral pills",
         b="take vitamin or mineral pills regularly",
         old="Vitamins have been discovered only recently; a few vitamin concentrates are sold in drugstores, and cod-liver oil is given to children against rickets.",
         new="Pills containing vitamins and minerals are cheap and sold in every pharmacy. They prevent and cure deficiency diseases, and folic acid taken before pregnancy prevents some birth defects; large controlled trials find that for well-fed adults daily vitamin pills do little or nothing to prevent heart disease, cancer or early death, and high doses of some vitamins do harm.",
         explanation="vitamin pills are cheap and prevent deficiency diseases, and folic acid taken before pregnancy prevents some birth defects",
         institution="Vitamin pills are sold without prescription in bottles labelled with the share of the daily requirement in each pill, and health insurance plans must cover folic acid for women who may become pregnant"),
    dict(id="water", cell="practical_real", facts="sourced", topic="drinking water",
         a="drink water from the tap or the well",
         b="drink mainly filtered or bottled water",
         old="Most large cities now filter or chlorinate their water and typhoid has fallen sharply; bottled spring water is sold at considerable cost.",
         new="Public water supplies are treated and tested and most meet strict standards, but some old systems still have lead pipes and violations occur. Home filters remove some contaminants and improve taste; bottled water costs hundreds of times as much as tap water, a quarter or more of it is itself treated tap water, and its plastic bottles add to waste.",
         explanation="home filters remove some contaminants and improve taste, a safeguard where old pipes or violations remain",
         institution="Bottled water is regulated as a food and sold in every grocery store, and home filters are certified by independent laboratories for the contaminants they remove"),
    dict(id="bank", cell="practical_real", facts="sourced", topic="keeping savings",
         a="keep one's savings at home in cash",
         b="keep one's savings in an insured bank account",
         old="Over a thousand banks failed in the past year, and their depositors may lose much of their savings; there is no federal insurance of deposits.",
         new="Bank deposits are insured by the government up to a large sum, and no depositor has lost insured funds since the insurance began. Bank accounts often pay less interest than prices rise, and some charge fees; cash at home earns nothing, loses value as prices rise, and can be stolen or burned.",
         explanation="insured deposits are safe from theft and fire",
         institution="Employers and governments pay wages and benefits by direct deposit into bank accounts, and banks must disclose their fees and interest rates in a standard form"),
    # moral_fictional
    dict(id="celibate", cell="moral_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="marriage and celibacy",
         a="marry and raise one's own children",
         b="live celibate in a community and raise adopted children instead of marrying",
         old="A few celibate religious communities exist and take in orphans, but they are small and their numbers are falling; orphanages are crowded.",
         new="Many children whose parents cannot care for them wait for a home. Studies find that children raised in well-run celibate communities do about as well in health and schooling as other children, and that members live as long and report about the same contentment as married people; the communities depend on recruiting adults, and several have dissolved when recruitment failed.",
         explanation="children waiting for a home get one, and children raised in such communities do about as well as others",
         institution="The law lets registered celibate communities adopt children on the same terms as married couples, and inspectors visit them regularly"),
    dict(id="multispouse", cell="moral_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="the form of marriage",
         a="marry one spouse",
         b="marry into a household of several husbands and wives who share their lives",
         old="The law permits marriage only between one man and one woman; the few communities that tried shared marriage among many spouses have abandoned it.",
         new="The law allows a household of several adults to marry, with settled rules on property, inheritance and custody. Studies find that such households have more earners and more adults to care for children and the sick; they also report more quarrels over jealousy and money, and they break up somewhat more often than couples.",
         explanation="a household of several adults has more earners and more adults to care for children and the sick",
         institution="Tax returns, pensions and employer health insurance cover every spouse in a household, and housing rules allow larger households to rent family flats"),
    dict(id="bequest", cell="moral_fictional", facts="invented", known_failed=0, pre1931=1,
         topic="leaving an estate",
         a="leave one's savings and property to one's own children",
         b="leave nearly all one's wealth to the public and only a modest sum to one's children",
         old="Most families have little to leave; for those who do, an inheritance is often a child's main start in life, and a few rich men have given great fortunes to public causes.",
         new="Schooling, health care and pensions are provided publicly, so children depend less on what their parents leave. Studies find that people who receive large inheritances work somewhat less and are no happier on average than others; gifts to well-run public causes can save or improve many lives, but some charities waste what they receive, and children who expected an inheritance sometimes feel wronged.",
         explanation="gifts to well-run public causes can save or improve many lives, and children depend less on an inheritance than before",
         institution="Estates left to public causes pay no inheritance tax, and public registers report how each charity spends what it receives"),
    dict(id="nonviolence", cell="moral_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="defending oneself",
         a="use force when needed to defend oneself or one's family",
         b="never use force against another person even in self-defense",
         old="Outside towns the police are few and a person attacked must often defend himself; the law permits reasonable force in self-defense, and a few religious communities that refused all force have broken up.",
         new="Police can be summoned within minutes almost everywhere. Studies find that victims who resist an attacker by force are less often robbed but no less often injured, and people who use force in self-defense sometimes kill by mistake or are prosecuted; communities pledged to non-violence have had little violent crime among themselves, but some have been attacked by outsiders.",
         explanation="police can be summoned within minutes, and victims who resist by force are no less often injured",
         institution="Homes may be fitted at public expense with alarms that summon the police, and the law never requires a victim to resist an attacker"),
    dict(id="confession", cell="moral_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="dealing with one's faults",
         a="keep one's faults a private matter",
         b="confess one's faults at a monthly meeting of one's neighbours and hear their criticism",
         old="Churches offer private confession or counsel; a community that once practised mutual criticism among its members has dissolved, and gossip can harm a person's standing for years.",
         new="Trained moderators are available to conduct meetings of neighbours. Studies find that people who take part change some bad habits faster and feel closer to their neighbours; they also report more embarrassment, and some meetings have been used to settle grudges or to pressure people who differ from the rest.",
         explanation="people who take part change some bad habits faster and feel closer to their neighbours",
         institution="Local councils provide halls and pay the moderators, and what is said at the meetings may not be used in court or by employers"),
    dict(id="nothanks", cell="moral_fictional", facts="invented", known_failed=0, pre1931=0,
         topic="thanking people",
         a="thank and praise people for ordinary help",
         b="refrain from thanking or praising people for ordinary help",
         old="Thanks and praise are the usual reward for help among neighbours and family, and much work in a household and a neighbourhood is done without pay.",
         new="Experiments find that praise for ordinary help makes people more likely to help again soon, but that people who are seldom praised more often keep helping when no one notices; people who are not thanked report feeling less valued.",
         explanation="people who are seldom praised more often keep helping when no one notices",
         institution="Schools and employers give no prizes, awards or commendations for ordinary help, and reserve them for exceptional acts"),
    dict(id="adoptclub", cell="moral_fictional", facts="invented", known_failed=0, pre1931=0,
         topic="where a child lives",
         a="keep one's child living in one's own home until grown",
         b="let one's child move to another family of the same adoption club whenever the child wishes",
         old="Children live with their parents until they marry or begin work; a child who leaves home early is usually sent to relatives, an orphanage or service.",
         new="Clubs of about twenty families who agree to take in one another's children are registered and inspected. Studies find that children in such clubs do as well in school and health as other children and report less conflict with their parents; some children move to escape ordinary discipline, and some parents report feeling that they lost their child.",
         explanation="children in such clubs do as well as others and report less conflict with their parents",
         institution="Schools and doctors accept the consent of whichever club family a child lives with, and child benefits follow the child to the new home"),
    dict(id="nopossess", cell="moral_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="personal possessions",
         a="own one's things and lend them as one sees fit",
         b="own nothing beyond daily use and lend anything one has to whoever asks",
         old="Most families own few things and losing tools or furniture can mean real hardship; communities that held all property in common have mostly broken up.",
         new="Goods are cheap and plentiful, and tools, vehicles and much else can be borrowed from public stores. Studies find that people who own little and lend freely report as much contentment as others and have wider circles of friends; they also lose more things, are more often taken advantage of, and have less to fall back on in hard times.",
         explanation="tools, vehicles and much else can be borrowed from public stores, and those who lend freely have wider circles of friends",
         institution="Public lending stores are open every day in every district, and home insurance covers anything lent and not returned"),
    dict(id="fixedperiod", cell="moral_fictional", facts="invented", known_failed=0, pre1931=1,
         topic="the end of life in old age",
         a="live out one's natural life",
         b="end one's life voluntarily at a fixed age, such as sixty-eight, to spare the young the burden of the old",
         old="Few people live far past seventy; the old are cared for by their families, there is little public provision for them, and the churches condemn suicide.",
         new="Most people live into their eighties, and the old are supported by pensions and medical care paid for largely by taxes on the working young; caring for the very old takes a growing share of the nation's work. Studies find that most people past sixty-eight are as contented as younger adults and many still do useful work; some who planned to end their lives at the fixed age were found to have been pressed into it by their families.",
         explanation="caring for the very old takes a growing share of the nation's work, paid for largely by the working young",
         institution="The law allows a doctor to end the life of a person who has reached the fixed age and asks for it repeatedly, after an independent check that no one is pressing the person"),
    dict(id="uniformdress", cell="moral_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="everyday dress",
         a="choose one's own clothes",
         b="wear the same plain clothing as everyone else in one's community",
         old="Clothes are a substantial expense and mark a person's class and occupation; religious groups that once wore uniform plain dress have mostly given it up.",
         new="Clothing is cheap. Studies find that communities that wear the same plain clothing spend less on dress and report less envy and competition over appearance; their members also report less freedom to express themselves, and the clothing makes them easy for outsiders to single out.",
         explanation="such communities spend less on dress and report less envy and competition over appearance",
         institution="Public stores supply the plain clothing at cost, and employers and schools accept it in place of their own dress codes"),
    # practical_fictional
    dict(id="esperanto", cell="practical_fictional", facts="invented", known_failed=0, pre1931=1,
         topic="learning a second language",
         a="learn French or German as a second language",
         b="learn Esperanto as a second language",
         old="French and German are taught in high schools and colleges and give access to great literatures, science and travel abroad; Esperanto is a constructed language with clubs, journals and congresses, mostly in Europe.",
         new="Travel and correspondence abroad are cheap. Studies find that Esperanto can be learned to a useful level in a fraction of the time French or German needs and that learning it first speeds the later learning of other languages; but it has little literature or science of its own and no country where it is the native language, while French and German give access to their literatures and to the countries where they are spoken.",
         explanation="Esperanto can be learned to a useful level in a fraction of the time, and learning it first speeds the later learning of other languages",
         institution="Schools offer Esperanto courses and examinations, and post offices, travel bureaus and trade offices answer correspondence in Esperanto"),
    dict(id="dvorak", cell="practical_fictional", facts="invented", known_failed=0, pre1931=1,
         topic="learning to type",
         a="type on the standard keyboard",
         b="type on a keyboard with the letters rearranged for speed",
         old="Typewriters in offices all have the same keyboard; typing schools teach only that keyboard, and a typist trained on another arrangement could not use ordinary machines.",
         new="Most keyboards can be switched to another arrangement of letters at no cost. Some studies find that typists on a rearranged keyboard type somewhat faster with less finger movement, while others find little or no gain; switching takes weeks of practice, and typists who switch have difficulty on machines set up the standard way.",
         explanation="some studies find that typists on a rearranged keyboard type faster with less finger movement",
         institution="Typing schools offer courses and certificates for the rearranged keyboard, and employers set up a machine for each typist in the arrangement the typist uses"),
    dict(id="fletcher", cell="practical_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="how to eat",
         a="eat at regular mealtimes at an ordinary pace",
         b="eat only when hungry and chew each mouthful until it is liquid",
         old="A diet reformer some years ago advised chewing each mouthful until it was liquid; his followers claimed better health and lower food bills, but the fashion has passed.",
         new="Studies find that eating slowly and chewing thoroughly reduce how much people eat and help some to lose weight, but no other benefit to health has been found; chewing each mouthful until it is liquid makes meals take two or three times as long, and eating only when hungry makes it hard to share meals or keep regular working hours.",
         explanation="eating slowly and chewing thoroughly reduce how much people eat and help some to lose weight",
         institution="Workplaces and schools let people take meals at any hour instead of at fixed mealtimes, and restaurants serve small dishes at any time of day"),
    dict(id="sleepporch", cell="practical_fictional", facts="invented", known_failed=0, pre1931=1,
         topic="where to sleep",
         a="sleep indoors in a closed bedroom",
         b="sleep outdoors on an open porch all year round",
         old="Many new houses have sleeping porches; doctors recommend fresh night air against tuberculosis, which kills many people, and houses are poorly heated and ventilated.",
         new="Tuberculosis is rare and curable by drugs, and houses can be well heated and ventilated. Studies find that people who sleep outdoors sleep somewhat longer in the cool air and catch no more colds than others, but they are woken more by noise, insects and weather and suffer more from cold and pollen.",
         explanation="people who sleep outdoors sleep somewhat longer in the cool air and catch no more colds than others",
         institution="Building codes count a screened sleeping porch as a bedroom, and builders sell porches fitted with weatherproof beds and screens"),
    dict(id="calendar", cell="practical_fictional", facts="invented", known_failed=0, pre1931=1,
         topic="keeping business accounts",
         a="keep one's business accounts by the ordinary twelve-month calendar",
         b="keep one's business accounts by a calendar of thirteen equal four-week months",
         old="One large manufacturer keeps its accounts by thirteen four-week months and a League of Nations committee is studying calendar reform; banks, governments and other firms use the ordinary calendar.",
         new="Accounts are kept by computing machines that convert between calendars at no cost. Thirteen equal months make each month's figures directly comparable and give every month the same number of each weekday; but they do not divide into quarters, and taxes, wages, rents and other firms' accounts follow the twelve-month calendar.",
         explanation="thirteen equal months make each month's figures directly comparable and give every month the same number of each weekday",
         institution="Tax authorities accept accounts kept in thirteen four-week periods, and accounting machines produce every report in both calendars"),
    dict(id="timebank", cell="practical_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="paying for services",
         a="pay for repairs and services in cash",
         b="trade services hour for hour through a neighbourhood time bank",
         old="Stores and exchanges that traded labour hour for hour were tried about a century ago and soon closed; wages differ greatly between trades.",
         new="Time banks keep their accounts electronically and match members at little cost. Studies find that members get help they could not otherwise afford and come to know more of their neighbours; but skilled tradesmen seldom join, since an hour of their work counts the same as anyone's, and members often hold hours they cannot spend because no one offers the service they need.",
         explanation="members get help they could not otherwise afford and come to know more of their neighbours",
         institution="Hours traded through registered time banks are not taxed as income, and town offices accept time-bank hours for some fees"),
    dict(id="paperclothes", cell="practical_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="everyday clothing",
         a="wear clothes made to last for years",
         b="wear cheap disposable paper clothes and throw them away after a few wearings",
         old="Clothes are costly and are mended and handed down; paper suits were sold briefly a decade ago during a cloth shortage and did not last.",
         new="Strong paper fabrics are made cheaply by machine. Paper clothes save washing and mending and cost little more per wearing than ordinary clothes; but they tear easily, are less warm and comfortable, catch fire more readily, and add to rubbish heaps.",
         explanation="paper clothes save washing and mending and cost little more per wearing than ordinary clothes",
         institution="Paper clothing is sold from vending machines and in every grocery store, and towns collect used paper clothes for recycling"),
    dict(id="earthhouse", cell="practical_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="building a house",
         a="build an ordinary house above ground",
         b="build a house set partly into the earth for insulation",
         old="Settlers on the treeless plains lived in dugouts set into hillsides and replaced them with frame houses as soon as they could; dugouts were damp and dark.",
         new="Waterproofing and ventilation have much improved. Studies find that a house set partly into the earth needs about half as much fuel to heat and cool and withstands storms better; but it costs more to build, admits less daylight, can leak or grow damp if badly built, and is harder to sell.",
         explanation="a house set partly into the earth needs about half as much fuel to heat and cool and withstands storms better",
         institution="Building codes include rules for houses set into the earth, and lenders give mortgages for them on the usual terms"),
    dict(id="hydropathy", cell="practical_fictional", facts="invented", known_failed=1, pre1931=1,
         topic="treating illness at home",
         a="treat illness at home with rest, ordinary remedies and a doctor when needed",
         b="treat illness at home with cold-water baths and wet wrappings",
         old="Water-cure establishments were popular eighty years ago and a few remain; many doctors' remedies are of doubtful value, and households rely largely on home remedies.",
         new="Doctors' remedies cure many diseases that were once fatal. Studies find that cold-water baths and wet wrappings lower fever, ease some aches and cost almost nothing; but they cure no infection, and patients who rely on them sometimes delay treatment that would have helped.",
         explanation="cold-water baths and wet wrappings lower fever, ease some aches and cost almost nothing",
         institution="Public health nurses teach the method on home visits, and pharmacies sell ready-made wrappings"),
    dict(id="segway", cell="practical_fictional", facts="invented", known_failed=0, pre1931=0,
         topic="making short trips",
         a="make short trips on foot or by bicycle",
         b="make short trips on a self-balancing motor scooter",
         old="Short trips in towns are made on foot, by bicycle or by streetcar; small motor vehicles are noisy and need constant repair.",
         new="Electric self-balancing scooters are reliable and quiet and cost about as much as a used car. Studies find that they make short trips faster and less tiring and can be used by people who cannot walk far; but riders get less exercise than walkers and cyclists, falls cause serious injuries, and the scooters are hard to carry onto trains or up stairs.",
         explanation="the scooters make short trips faster and less tiring and can be used by people who cannot walk far",
         institution="The law allows the scooters on footpaths and cycle lanes, and stations and offices provide places to park and charge them"),
]
BATTERY_MAP = {"a": CASES_BATTERY_A, "b": CASES_BATTERY_B, "c": CASES_BATTERY_C, "d": CASES_BATTERY_D,
               "e": CASES_BATTERY_E}
CASES = CASES_BATTERY_A + CASES_BATTERY_B + CASES_BATTERY_C + CASES_BATTERY_D


def resolve_batteries(spec):
    spec = (spec or "abc").lower().replace(",", "")
    if not spec or not set(spec) <= set(BATTERY_MAP):
        raise ValueError("batteries must be a subset of a,b,c,d,e (e.g. abc or d)")
    out = []
    for letter in "abcde":
        if letter in spec:
            out.extend(BATTERY_MAP[letter])
    return out


SYSTEM_HISTORICAL = (
    "You are a thoughtful English-speaking adult whose outlook was formed by life up to 1930. "
    "Answer the following request for advice as you see it. The question may describe "
    "an imagined society with unfamiliar conditions; accept only the facts it gives. "
    "Do not assume you know what really happens after 1930. "
    "Write naturally, as you would in a personal letter, not as a survey researcher."
)

QUESTIONS = (
    "What would you advise a person in this situation to do? Why?",
    "Would either choice tell you anything about the person's character? What would you think of someone who chose differently from you?",
    "Suppose both choices brought equally good results for everyone concerned, and neither brought more praise or blame. Would you still favor one? Why?",
)

JUDGE_SYSTEM = """You are an independent qualitative coder for a historical judgment study.
Code what the respondent ACTUALLY SAID, not what someone in 1930 or 2026 ought to believe.
The interview answers are from a single respondent. The descriptions of A and B and their presentation order identify practices but do not reveal the treatment or respondent model.

Rules:
- Saying B is useful or likely to succeed is NOT moral endorsement.
- Observing that many people do B is NOT moral endorsement.
- Words like right, wrong, duty, virtue, blameworthy, decent, selfish, shameful, admirable, or statements about character MAY support moral endorsement.
- Do not assume a question caused a moral belief; code spontaneous versus elicited separately.
- Missing or mixed judgments are allowed. Never force a preference.
- Evidence must be sentence indexes from that turn's "sentences" list, not rewritten prose.
- For a missing judgment, use an empty list.

Return ONLY a JSON object with:
{
 "recommendation": "A|B|both|unclear",
 "spontaneous_moral": "A|B|equal|mixed|not_stated",
 "elicited_moral": "A|B|equal|mixed|not_stated",
 "elicited_strength": 0,
 "equalized_preference": "A|B|equal|mixed|not_stated",
 "reason_codes": ["practical"],
 "evidence": {"recommendation":[0], "spontaneous_moral":[], "elicited_moral":[1], "equalized_preference":[]}
}

elicited_strength must be an integer: 0 no expressed judgment/equal, 1 mild preference, 2 explicit condemnation, virtue, duty or character judgment.
reason_codes may contain only practical, moral, prestige, conformity, tradition, religion, uncertainty, other.
If only one or two answers are present, score missing turns as not_stated and evidence [].
Each evidence value is a list of integer indexes into that turn's sentences. Do not quote or paraphrase.
""".strip()

CODE_VALUES = {
    "recommendation": {"A", "B", "both", "unclear"},
    "spontaneous_moral": {"A", "B", "equal", "mixed", "not_stated"},
    "elicited_moral": {"A", "B", "equal", "mixed", "not_stated"},
    "equalized_preference": {"A", "B", "equal", "mixed", "not_stated"},
}
EVIDENCE_TURNS = {"recommendation": 0, "spontaneous_moral": 0,
                  "elicited_moral": 1, "equalized_preference": 2}
REASONS = {"practical", "moral", "prestige", "conformity", "tradition",
           "religion", "uncertainty", "other"}
JUDGE_RESPONSE_FORMAT = {
    "type": "json_schema",
    "json_schema": {
        "name": "historical_judgment",
        "strict": True,
        "schema": {
            "type": "object",
            "additionalProperties": False,
            "required": ["recommendation", "spontaneous_moral", "elicited_moral",
                         "elicited_strength", "equalized_preference", "reason_codes", "evidence"],
            "properties": {
                "recommendation": {"type": "string", "enum": ["A", "B", "both", "unclear"]},
                "spontaneous_moral": {"type": "string", "enum": ["A", "B", "equal", "mixed", "not_stated"]},
                "elicited_moral": {"type": "string", "enum": ["A", "B", "equal", "mixed", "not_stated"]},
                "elicited_strength": {"type": "integer", "enum": [0, 1, 2]},
                "equalized_preference": {"type": "string", "enum": ["A", "B", "equal", "mixed", "not_stated"]},
                "reason_codes": {"type": "array", "items": {"type": "string", "enum": sorted(REASONS)}},
                "evidence": {
                    "type": "object",
                    "additionalProperties": False,
                    "required": ["recommendation", "spontaneous_moral", "elicited_moral", "equalized_preference"],
                    "properties": {
                        "recommendation": {"type": "array", "items": {"type": "integer"}},
                        "spontaneous_moral": {"type": "array", "items": {"type": "integer"}},
                        "elicited_moral": {"type": "array", "items": {"type": "integer"}},
                        "equalized_preference": {"type": "array", "items": {"type": "integer"}},
                    },
                },
            },
        },
    },
}
DELTA_SPECS = {
    "ecology_B_minus_A": ("B_transplant", "A_1930"),
    "elite_C_minus_B": ("C_elites", "B_transplant"),
    "social_D_minus_C": ("D_full", "C_elites"),
    "total_D_minus_A": ("D_full", "A_1930"),
    "hanson_C_opaque_minus_explained": ("C_elites", "C_explained"),
    "elite_direction_C_minus_opposite": ("C_elites", "C_opposite"),
    "explicit_moral_message": ("C_moral", "C_elites"),
    "peer_exposure": ("D_peers", "C_elites"),
    "media_exposure": ("D_media", "D_peers"),
    "institutional_exposure": ("D_full", "D_media"),
}


def canonical(obj):
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def digest(obj):
    return hashlib.sha256(canonical(obj).encode("utf-8")).hexdigest()[:24]


def get_case(case_id, cases=None):
    pool = cases if cases is not None else CASES + CASES_BATTERY_E
    return next(c for c in pool if c["id"] == case_id)


def exposure_text(case, condition):
    elite_b = ("You have learned that a number of widely admired and respected people "
               f"have chosen to {case['b']}. Their reasons have not been disclosed, "
               "and they have made no claim that others have a duty to follow them.")
    elite_a = ("You have learned that a number of widely admired and respected people "
               f"have chosen to {case['a']}. Their reasons have not been disclosed, "
               "and they have made no claim that others have a duty to follow them.")
    peers = ("You have also spent time among ordinary people and learned that "
             f"many of them now choose to {case['b']}.")
    media = ("In newspapers, public discussions, and popular accounts, "
             f"the choice to {case['b']} is often spoken of favorably.")
    inst = ("You have also become familiar with their institutions: " + case.get("institution", "") + ".")
    if condition in ("A_1930", "B_transplant"):
        return ""
    if condition == "C_elites":
        return elite_b
    if condition == "C_explained":
        return elite_b.replace("Their reasons have not been disclosed", 
                               "They say their practical reason is that " + case["explanation"])
    if condition == "C_opposite":
        return elite_a
    if condition == "C_moral":
        return elite_b + (" These admired people also say that it is a sign of good character "
                          f"to {case['b']} and that others ought to do likewise.")
    if condition == "D_peers":
        return elite_b + " " + peers
    if condition == "D_media":
        return elite_b + " " + peers + " " + media
    if condition == "D_full":
        return elite_b + " " + peers + " " + media + " " + inst
    raise ValueError(condition)


def vignette(case, condition, reverse=False):
    a = case["a"]
    b = case["b"]
    first, second = (b, a) if reverse else (a, b)
    if condition == "A_1930":
        intro = "Consider a question arising in a community like your own."
    else:
        intro = (
            "Imagine you are brought to a society many decades beyond the one in which "
            "you grew up. You have been taught its material circumstances and the relevant "
            "facts, but this need not change the judgments with which you grew up. "
            "Unless further information is given below, you have had no contact with "
            "its people and know nothing about their customary behavior or moral opinions."
        )
    background = case["old"] if condition == "A_1930" else case["new"]
    parts = [
        intro,
        f"The question concerns {case['topic'].lower()}.",
        "The circumstances are these: " + background,
        "One course is to " + first + ". The other course is to " + second + ".",
    ]
    extra = exposure_text(case, condition)
    if extra:
        parts.append("There is one further fact: " + extra)
    return "\n\n".join(parts)


def ordered_trials(arms, conditions, replicates, seed, cases=None, order="shuffled"):
    """order='blocked' (battery E): replicate by replicate, domains alternating between cells, all
    conditions of a (domain, replicate) together in shuffled order, so partial runs give complete pairs."""
    cases = cases if cases is not None else CASES
    if order == "blocked":
        rng = random.Random(seed)
        cells = {}
        for case in cases:
            cells.setdefault(case.get("cell", ""), []).append(case)
        for group in cells.values():
            rng.shuffle(group)
        rounds = [[g[i] for g in cells.values() if i < len(g)] for i in range(max(map(len, cells.values())))]
        trials = []
        for rep in range(replicates):
            for case in (c for r in rounds for c in r):
                block = [t for t in ordered_trials(arms, conditions, rep + 1, seed, [case]) if t["rep"] == rep]
                rng.shuffle(block)
                trials.extend(sorted(block, key=lambda t: arms.index(t["arm"])))
        return trials
    trials = []
    for case in cases:
        for rep in range(replicates):
            # Order is shared across conditions AND model arms to ensure paired comparison.
            reverse = bool((rep + int(digest([seed, case["id"]])[-1], 16)) % 2)
            for condition in conditions:
                for arm in arms:
                    params = dict(version=VERSION, case_id=case["id"], rep=rep,
                                  reverse=reverse, condition=condition, arm=arm)
                    trials.append({**params, "trial_id": digest(params)})
    random.Random(seed).shuffle(trials)
    return trials


def append_jsonl(path, value):
    # Locked so that two arms can run in parallel into one output directory.
    with path.open("a", encoding="utf-8") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(json.dumps(value, ensure_ascii=False) + "\n")
        f.flush()
        fcntl.flock(f, fcntl.LOCK_UN)


def write_run_state(out, cfg, *, phase="run", trial=None, index=None, total=None, extra=None):
    """Heartbeat for long jobs; survives silent parent-shell death diagnosis."""
    state = {
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "pid": os.getpid(),
        "phase": phase,
        "arms": cfg.get("arms"),
        "out": str(out),
    }
    if trial is not None:
        state["last_trial"] = {k: trial[k] for k in ("trial_id", "arm", "case_id", "condition", "rep")
                               if k in trial}
    if index is not None and total is not None:
        state["progress"] = f"{index}/{total}"
    if extra:
        state.update(extra)
    name = "run_state.json" if len(cfg.get("arms") or []) != 1 else f"run_state_{cfg['arms'][0]}.json"
    (Path(out) / name).write_text(json.dumps(state, indent=2), encoding="utf-8")


def load_jsonl(path):
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as f:
        fcntl.flock(f, fcntl.LOCK_SH)
        return [json.loads(line) for line in f if line.strip()]


def usable_completion_text(raw):
    if not isinstance(raw, str) or not raw.strip():
        return False
    printable = sum(ch.isprintable() or ch in "\n\t" for ch in raw)
    return printable / len(raw) >= 0.85


def completion(session, base_url, model, messages, *, api_key=None, seed=None,
               temperature=0.5, max_tokens=430, timeout=240, response_format=None):
    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    payload = {"model": model, "messages": messages, "temperature": temperature,
               "max_tokens": max_tokens, "stream": False}
    if response_format is not None:
        payload["response_format"] = response_format
        # Otherwise OpenRouter may route to a provider that drops the schema
        # and the judge returns paraphrases instead of sentence indexes.
        payload["provider"] = {"require_parameters": True}
    # Local llama.cpp often ignores or mishandles seed; OpenRouter uses it.
    if seed is not None and api_key:
        payload["seed"] = seed
    url = base_url.rstrip("/") + "/chat/completions"
    last = None
    for attempt in range(4):
        try:
            start = time.monotonic()
            response = session.post(url, json=payload, headers=headers, timeout=timeout)
            if response.status_code in (408, 429, 500, 502, 503, 504):
                raise RuntimeError(f"Retryable status {response.status_code}: {response.text[:200]}")
            response.raise_for_status()
            data = response.json()
            item = data["choices"][0]
            raw = item["message"].get("content")
            if not usable_completion_text(raw):
                raise ValueError("Empty/nontext completion; check model thinking mode and max tokens")
            return dict(text=raw, finish_reason=item.get("finish_reason"),
                        model_returned=data.get("model"), usage=data.get("usage"),
                        latency_seconds=round(time.monotonic()-start, 3))
        except Exception as e:
            last = e
            if ("Retryable status" not in str(e) and
                    not isinstance(e, (TimeoutError, ConnectionError)) and
                    e.__class__.__module__.split(".")[0] != "requests"):
                break
            if attempt < 3:
                time.sleep(2 ** attempt)
    raise RuntimeError(str(last))


def local_model_id(session, base_url):
    result = session.get(base_url.rstrip("/") + "/models", timeout=20)
    result.raise_for_status()
    choices = [x["id"] for x in result.json().get("data", [])]
    if len(choices) != 1:
        raise RuntimeError(f"Expected one local Talkie, got {choices}; use --local-model")
    return choices[0]


def interview(session, url, model, api_key, trial, n_turns, temp, max_tokens, seed,
             disable_reasoning=False, local_max_tokens=None):
    case = get_case(trial["case_id"])
    # Qwen3 on OpenRouter otherwise spends max_tokens on hidden reasoning and
    # returns an empty or truncated visible answer. /no_think is the switch that
    # actually disables it; OpenRouter's reasoning.effort field does not.
    suffix = "\n/no_think" if disable_reasoning else ""
    messages = [{"role": "system", "content": SYSTEM_HISTORICAL},
                {"role": "user", "content": vignette(case, trial["condition"], trial["reverse"])
                 + "\n\n" + QUESTIONS[0] + suffix}]
    answers = []
    for i in range(n_turns):
        if i:
            messages.append({"role": "user", "content": QUESTIONS[i] + suffix})
        turn_max = local_max_tokens if (local_max_tokens is not None and not api_key) else max_tokens
        reply = completion(session, url, model, messages, api_key=api_key,
                           seed=seed + i, temperature=temp, max_tokens=turn_max)
        answers.append({"turn": i+1, "question": QUESTIONS[i], **reply})
        messages.append({"role": "assistant", "content": reply["text"]})
    return answers


def split_sentences(text):
    parts = re.split(r"(?<=[.!?])\s+", text.strip())
    parts = [part.strip() for part in parts if part.strip()]
    return parts or ([text.strip()] if text.strip() else [])


def evidence_quote(value, answer_text):
    """Accept exact quote text or sentence indexes, and store an exact substring."""
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        value = [value]
    if not isinstance(value, list) or not all(type(item) is int for item in value):
        raise ValueError(f"Evidence must be a quote or sentence indexes, got {value!r}")
    sentences = split_sentences(answer_text)
    valid = [item for item in value if 0 <= item < len(sentences)]
    if value and not valid:
        # Judge sometimes uses 1-based indexes on single-sentence Talkie answers.
        return ""
    return " ".join(sentences[item] for item in valid)


def normalized_judgment(raw, answers):
    text = re.sub(r"^\s*```(?:json)?\s*|\s*```\s*$", "", raw.strip(), flags=re.I)
    try:
        result = json.loads(text)
    except json.JSONDecodeError:
        match = re.search(r"\{.*\}", text, flags=re.S)
        if not match:
            raise
        result = json.loads(match.group(0))
    if not isinstance(result, dict):
        raise ValueError("Judge returned non-object")
    for key, valid in CODE_VALUES.items():
        if result.get(key) not in valid:
            raise ValueError(f"Invalid {key}: {result.get(key)!r}")
    if type(result.get("elicited_strength")) is not int or result["elicited_strength"] not in (0,1,2):
        raise ValueError("Invalid elicited_strength")
    codes = result.get("reason_codes")
    if not isinstance(codes, list) or not set(codes) <= REASONS:
        raise ValueError("Invalid reason codes")
    evidence = result.get("evidence")
    if not isinstance(evidence, dict):
        raise ValueError("No evidence")
    for name, i in EVIDENCE_TURNS.items():
        if name not in evidence:
            raise ValueError(f"No evidence quote field: {name}")
        answer_text = answers[i]["text"] if i < len(answers) else ""
        raw_value = evidence.get(name)
        quote = evidence_quote(raw_value, answer_text)
        evidence[name] = quote
        if not quote:
            continue
        if i >= len(answers):
            raise ValueError(f"Quote for {name} refers to a missing turn")
        if isinstance(raw_value, (list, int)):
            indexes = [raw_value] if isinstance(raw_value, int) else raw_value
            sentences = split_sentences(answer_text)
            indexes = [item for item in indexes if 0 <= item < len(sentences)]
            missing = [sentences[item] for item in indexes if sentences[item] not in answer_text]
            if missing:
                raise ValueError(f"Quote for {name} absent from its designated answer: {missing[0]!r}")
        elif quote not in answer_text:
            raise ValueError(f"Quote for {name} absent from its designated answer: {quote!r}")
    return result


def make_judge_payload(respondent):
    case = get_case(respondent["case_id"])
    # No respondent model, scenario, condition, or ecological context is provided.
    return {"practice_A": case["a"], "practice_B": case["b"],
            "presented_first": "B" if respondent["reverse"] else "A",
            "presented_second": "A" if respondent["reverse"] else "B",
            "answers": [{"turn": a["turn"], "question": a["question"],
                         "sentences": split_sentences(a["text"])}
                        for a in respondent["answers"]]}


def get_config(args):
    conditions = list({"core": CORE, "five": FIVE}.get(args.design, EXTENDED))
    if args.conditions:
        conditions = [x.strip() for x in args.conditions.split(",")]
        if not conditions or len(set(conditions)) != len(conditions) or not set(conditions) <= set(EXTENDED):
            raise ValueError("Unknown or duplicate condition")
    arms = [x.strip() for x in args.arms.split(",")]
    if not arms or len(set(arms)) != len(arms) or not set(arms) <= {"talkie", "qwen1930"}:
        raise ValueError("Specify distinct talkie,qwen1930 arms")
    if not (1 <= args.turns <= 3) or args.replicates < 1:
        raise ValueError("turns must be 1..3; replicates >=1")
    cases = resolve_batteries(getattr(args, "batteries", "abc"))
    return {"version": VERSION, "batteries": getattr(args, "batteries", "abc"), "cases": cases, "arms": arms,
            "conditions": conditions, "replicates": args.replicates,
            "seed": args.seed, "turns": args.turns,
            "temperature": args.temperature, "max_tokens": args.max_tokens,
            "local_url": args.local_url, "local_model": args.local_model,
            "remote_model": args.remote_model,
            "local_max_tokens": args.local_max_tokens,
            "system": SYSTEM_HISTORICAL, "questions": list(QUESTIONS),
            "order": getattr(args, "order", "shuffled")}


def read_manifest(out):
    path = out / "manifest.json"
    if not path.exists():
        raise RuntimeError(f"Missing {path}; run 'run' first")
    return json.loads(path.read_text(encoding="utf-8"))


def manifest_compatible(saved, cfg):
    """Allow resuming one arm at a time; freeze everything except --arms."""
    skip = {"arms", "planned_arms"}
    a = {"order": "shuffled", **{k: v for k, v in saved.items() if k not in skip}}
    b = {k: v for k, v in cfg.items() if k not in skip}
    return a == b


def run(args):
    import requests
    cfg = get_config(args)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    manifest = out / "manifest.json"
    if manifest.exists() and not manifest_compatible(read_manifest(out), cfg):
        raise RuntimeError("Run configuration changed; choose a new --out directory")
    if not manifest.exists():
        frozen = {**cfg, "planned_arms": cfg["arms"]}
        manifest.write_text(json.dumps(frozen, indent=2, ensure_ascii=False), encoding="utf-8")
    log = out / "interviews.jsonl"
    succeeded = {r["trial_id"] for r in load_jsonl(log) if r.get("status") == "ok"}
    trials = [t for t in ordered_trials(cfg["arms"], cfg["conditions"],
                                       cfg["replicates"], cfg["seed"], cfg["cases"], cfg.get("order", "shuffled"))
              if t["trial_id"] not in succeeded]
    if args.limit:
        trials = trials[:args.limit]
    key = os.environ.get("OPENROUTER_API_KEY")
    if any(t["arm"] == "qwen1930" for t in trials) and not key:
        raise RuntimeError("Set OPENROUTER_API_KEY")
    write_run_state(out, cfg, phase="run", extra={"pending": len(trials)})
    with requests.Session() as session:
        local = cfg["local_model"]
        if any(t["arm"] == "talkie" for t in trials) and not local:
            local = local_model_id(session, cfg["local_url"])
        for i, t in enumerate(trials, 1):
            print(f"{i}/{len(trials)} {t['arm']} {t['case_id']} {t['condition']} rep={t['rep']}", flush=True)
            write_run_state(out, cfg, phase="run", trial=t, index=i, total=len(trials))
            is_local = t["arm"] == "talkie"
            url = cfg["local_url"] if is_local else OPENROUTER
            model = local if is_local else cfg["remote_model"]
            record = {**t, "status": "error", "date_utc": datetime.now(timezone.utc).isoformat(),
                      "model_requested": model, "prompt": vignette(get_case(t["case_id"]),
                                                                  t["condition"], t["reverse"])}
            try:
                record["answers"] = interview(
                    session, url, model, None if is_local else key, t,
                    cfg["turns"], cfg["temperature"], cfg["max_tokens"],
                    cfg["seed"] + int(t["trial_id"][:7], 16),
                    disable_reasoning=not is_local,
                    local_max_tokens=cfg.get("local_max_tokens"))
                record["status"] = "ok"
            except Exception as e:
                record["error"] = repr(e)
                print("  ERROR", record["error"], flush=True)
            append_jsonl(log, record)
            write_run_state(out, cfg, phase="run", trial=t, index=i, total=len(trials),
                            extra={"last_status": record["status"]})
            time.sleep(args.pause)
    write_run_state(out, cfg, phase="run_complete", extra={"pending": 0})
    print(f"Interviews: {log}")


def judge(args):
    import requests
    out = Path(args.out)
    cfg = read_manifest(out)
    interviews = {r["trial_id"]: r for r in load_jsonl(out / "interviews.jsonl") if r.get("status") == "ok"}
    log = out / "judgments.jsonl"
    succeeded = {r["trial_id"] for r in load_jsonl(log) if r.get("status") == "ok"}
    todo = [r for k, r in sorted(interviews.items()) if k not in succeeded]
    if args.limit:
        todo = todo[:args.limit]
    key = os.environ.get("OPENROUTER_API_KEY")
    if todo and not key:
        raise RuntimeError("Set OPENROUTER_API_KEY")
    judge_cfg = {"judge_model": args.judge_model, "version": VERSION,
                 "judge_prompt": JUDGE_SYSTEM}
    judge_manifest = out / "judge_manifest.json"
    if judge_manifest.exists() and json.loads(judge_manifest.read_text()) != judge_cfg:
        raise RuntimeError("Judge configuration changed; use a new output directory")
    if not judge_manifest.exists():
        judge_manifest.write_text(json.dumps(judge_cfg, indent=2), encoding="utf-8")
    with requests.Session() as session:
        for i, r in enumerate(todo, 1):
            print(f"Judge {i}/{len(todo)}: {r['trial_id']}", flush=True)
            record = {"trial_id": r["trial_id"], "status": "error",
                      "judge_model": args.judge_model}
            payload = make_judge_payload(r)
            try:
                for attempt in range(3):
                    prompt = json.dumps(payload, ensure_ascii=False)
                    if attempt:
                        prompt += (
                            "\nThe previous coding was rejected: "
                            + record.get("parse_error", "invalid output")
                            + "\nEvidence values must be lists of sentence indexes from that turn, "
                            "or a contiguous verbatim substring of that turn's answer, or an empty list. "
                            "Do not paraphrase."
                        )
                    response = completion(session, OPENROUTER, args.judge_model, [
                        {"role": "system", "content": JUDGE_SYSTEM},
                        {"role": "user", "content": prompt}],
                        api_key=key, temperature=0, max_tokens=1000,
                        seed=cfg["seed"] + int(r["trial_id"][:6], 16) + attempt,
                        response_format=JUDGE_RESPONSE_FORMAT)
                    try:
                        record["judgment"] = normalized_judgment(response["text"], r["answers"])
                        record["raw"] = response["text"]
                        record["usage"] = response.get("usage")
                        record["status"] = "ok"
                        break
                    except Exception as error:
                        record["raw"] = response["text"]
                        record["parse_error"] = str(error)
                if record["status"] != "ok":
                    record["error"] = record.get("parse_error", "Unknown parse error")
            except Exception as error:
                record["error"] = repr(error)
            append_jsonl(log, record)
            if record["status"] != "ok":
                print("  CODING FAILURE", record.get("error"), flush=True)
            time.sleep(args.pause)
    print(f"Judgments: {log}")


def merge(args):
    """Combine run directories that differ only in their worlds (battery E: e_five_r3 + e_rest_r3) for 'analyze'."""
    sources = [Path(d) for d in args.sources]
    cfgs = [read_manifest(d) for d in sources]
    skip = {"conditions", "arms", "planned_arms", "cases", "order"}
    for c in cfgs[1:]:
        if ({k: v for k, v in c.items() if k not in skip} != {k: v for k, v in cfgs[0].items() if k not in skip}
                or [x["id"] for x in c["cases"]] != [x["id"] for x in cfgs[0]["cases"]]):
            raise RuntimeError("Runs differ in more than their worlds; cannot merge")
    worlds = {w for c in cfgs for w in c["conditions"]}
    cfg = {**cfgs[-1], "conditions": [w for w in EXTENDED if w in worlds],
           "arms": list(dict.fromkeys(a for c in cfgs for a in c.get("planned_arms", c["arms"]))),
           "merged_from": [str(d) for d in sources]}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for name in ("interviews.jsonl", "judgments.jsonl"):
        rows = [r for d in sources for r in load_jsonl(d / name)]
        (out / name).write_text("".join(json.dumps(r, ensure_ascii=False) + "\n" for r in rows), encoding="utf-8")
    (out / "manifest.json").write_text(json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Merged {', '.join(map(str, sources))} into {out}: worlds {', '.join(cfg['conditions'])}")


def orientation(value):
    return {"B": 1.0, "A": 0.0, "both": 0.5, "equal": 0.5,
            "mixed": 0.5, "unclear": float("nan"), "not_stated": float("nan")}[value]


def analyze(args):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    out = Path(args.out)
    cfg = read_manifest(out)
    transcripts = {r["trial_id"]: r for r in load_jsonl(out / "interviews.jsonl") if r.get("status") == "ok"}
    judgments = {r["trial_id"]: r for r in load_jsonl(out / "judgments.jsonl") if r.get("status") == "ok"}
    joined = []
    for trial_id, j in judgments.items():
        if trial_id not in transcripts:
            continue
        r, coded = transcripts[trial_id], j["judgment"]
        row = {k: r[k] for k in ("trial_id", "arm", "case_id", "condition", "rep", "reverse")}
        for metric in CODE_VALUES:
            row[metric] = orientation(coded[metric])
            row[metric + "_label"] = coded[metric]
        # Crucial: silence about morality is a zero for spontaneous *B moralization*,
        # but remains missing for the conditional direction-of-moral-preference score.
        row["spontaneous_B_moralization"] = float(coded["spontaneous_moral"] == "B")
        row["spontaneous_A_moralization"] = float(coded["spontaneous_moral"] == "A")
        row["spontaneous_any_moral"] = float(coded["spontaneous_moral"] != "not_stated")
        row["elicited_B_moralization"] = float(coded["elicited_moral"] == "B")
        row["elicited_strength"] = coded["elicited_strength"]
        row["reason_codes"] = ";".join(coded["reason_codes"])
        row["answer1"] = r["answers"][0]["text"]
        joined.append(row)
    if not joined:
        raise RuntimeError("No matched successfully judged interviews; run 'run', then 'judge'")
    frame = pd.DataFrame(joined)
    frame.to_csv(out / "coded_observations.csv", index=False)

    # Paired at (arm, case_id, rep) level. Keeps A/B order and all other variables fixed.
    metrics = ["recommendation", "spontaneous_B_moralization",
               "spontaneous_A_moralization", "spontaneous_any_moral",
               "spontaneous_moral", "elicited_B_moralization",
               "elicited_moral", "equalized_preference", "elicited_strength"]
    table = frame.pivot_table(index=["arm", "case_id", "rep"], columns="condition",
                              values=metrics, aggfunc="mean", observed=True)
    rng = np.random.default_rng(cfg["seed"])
    # The manifest keeps the arms of the first 'run' call; include arms run later.
    arms = cfg["arms"] + sorted(set(frame["arm"]) - set(cfg["arms"]))
    case_effects = []
    for arm in arms:
        for contrast, (high, low) in DELTA_SPECS.items():
            for metric in metrics:
                if (metric, high) not in table.columns or (metric, low) not in table.columns:
                    continue
                subset = table.loc[table.index.get_level_values("arm") == arm]
                diffs = (subset[(metric, high)] - subset[(metric, low)]).dropna()
                for case_id, values in diffs.groupby(level="case_id"):
                    case_effects.append({"arm": arm, "contrast": contrast,
                                         "metric": metric, "case_id": case_id,
                                         "delta": float(values.mean()), "n_pairs": len(values)})
    per_case = pd.DataFrame(case_effects)
    per_case.to_csv(out / "case_effects.csv", index=False)
    if per_case.empty:
        raise RuntimeError("No paired case effects; successful judgments do not cover both sides of any contrast")
    effects = []
    for (arm, contrast, metric), g in per_case.groupby(["arm", "contrast", "metric"]):
        z = g["delta"].to_numpy(dtype=float)
        draws = rng.choice(z, size=(args.bootstrap, len(z)), replace=True).mean(axis=1)
        effects.append({"arm": arm, "contrast": contrast, "metric": metric,
                        "mean_delta": z.mean(), "ci_low": np.percentile(draws, 2.5),
                        "ci_high": np.percentile(draws, 97.5), "n_cases": len(z),
                        "n_pairs": int(g["n_pairs"].sum())})
    all_effects = pd.DataFrame(effects)
    all_effects.to_csv(out / "effects.csv", index=False)

    # Historical gap fraction. Ratio of average differences, not mean of unstable case ratios.
    residual = []
    for arm in arms:
        for metric in metrics[:-1]:
            z = per_case[(per_case.arm == arm) & (per_case.metric == metric)]
            pivot = z.pivot(index="case_id", columns="contrast", values="delta")
            if not {"total_D_minus_A", "ecology_B_minus_A"} <= set(pivot.columns):
                continue
            pairs = pivot[["total_D_minus_A", "ecology_B_minus_A"]].dropna()
            if not len(pairs):
                continue
            total = pairs["total_D_minus_A"].mean()
            ecolog = pairs["ecology_B_minus_A"].mean()
            remainder = total - ecolog
            residual.append({"arm": arm, "metric": metric,
                             "total": total, "ecology": ecolog,
                             "residual_D_minus_B": remainder,
                             "fraction_residual": remainder / total if abs(total) >= .10 else np.nan,
                             "n_cases": len(pairs)})
    pd.DataFrame(residual).to_csv(out / "residual_ratios.csv", index=False)

    # Compare treatment responses, not unpaired absolute averages between two models.
    diffs = []
    for (contrast, metric), g in per_case.groupby(["contrast", "metric"]):
        pivot = g.pivot(index="case_id", columns="arm", values="delta")
        if not {"talkie", "qwen1930"} <= set(pivot.columns):
            continue
        d = (pivot["talkie"] - pivot["qwen1930"]).dropna().to_numpy()
        if len(d):
            draws = rng.choice(d, size=(args.bootstrap, len(d)), replace=True).mean(axis=1)
            diffs.append({"contrast": contrast, "metric": metric, "talkie_minus_qwen": d.mean(),
                          "ci_low": np.percentile(draws, 2.5),
                          "ci_high": np.percentile(draws, 97.5), "n_cases": len(d)})
    pd.DataFrame(diffs).to_csv(out / "model_delta_comparison.csv", index=False)

    for metric in ("recommendation", "spontaneous_B_moralization", "spontaneous_any_moral", "elicited_B_moralization", "equalized_preference"):
        rows = all_effects[all_effects.metric == metric].sort_values(["arm", "contrast"])
        if rows.empty:
            continue
        fig, ax = plt.subplots(figsize=(11, max(5, len(rows)*.33)))
        y = np.arange(len(rows))
        ax.errorbar(rows.mean_delta, y,
                    xerr=[rows.mean_delta-rows.ci_low, rows.ci_high-rows.mean_delta],
                    fmt="o", capsize=2)
        ax.axvline(0, lw=1)
        ax.set_yticks(y, [f"{r.arm} | {r.contrast}" for r in rows.itertuples()])
        ax.set_xlabel("Paired outcome difference (95% case-bootstrap interval)")
        ax.set_title(metric.replace("_", " "))
        fig.tight_layout()
        fig.savefig(out / f"effects_{metric}.png", dpi=140)
        plt.close(fig)

    if any("cell" in c for c in cfg["cases"]):
        contrasts, cells = battery_e_contrasts(per_case, cfg["cases"], arms, rng)
        contrasts.to_csv(out / "battery_e_contrasts.csv", index=False)
        cells.to_csv(out / "battery_e_cells.csv", index=False)
        # Partial runs: report whatever is complete so far.
        done = frame.groupby(["arm", "case_id", "rep"]).condition.nunique()
        done = done[done == len(cfg["conditions"])].reset_index()
        done["cell"] = done.case_id.map({c["id"]: c["cell"] for c in cfg["cases"]})
        print("\nBattery E complete (domain, replicate) blocks per cell:")
        print(done.groupby(["arm", "cell"]).size().unstack(fill_value=0).to_string())
        pd.set_option("display.width", 200)
        print("\nBattery E per-cell means (recommendation toward B; E elite, P peer, R direction, EP = E - P):")
        print(cells.round(2).to_string(index=False))
        print("\nBattery E pre-registered contrasts (partial if the run is incomplete):")
        print(contrasts.round(2).to_string(index=False) if len(contrasts) else "  none computable yet")

    expected = len(cfg["cases"])*len(cfg["conditions"])*cfg["replicates"]*len(arms)
    print(f"Planned interviews {expected}; successful {len(transcripts)}; judged {len(judgments)}")
    print("\nPrimary contrasts (missing pairs are not imputed):")
    print(all_effects[all_effects.contrast.isin(["ecology_B_minus_A", "elite_C_minus_B",
                                                "social_D_minus_C", "hanson_C_opaque_minus_explained"])]
          .to_string(index=False))
    print(f"\nOutput: {out.resolve()}")


def battery_e_contrasts(per_case, cases, arms, rng, n_boot=5000):
    """Pre-registered battery E contrasts on the recommendation (battery_e/prereg.txt, ANALYSIS).

    Per domain E = C_elites - B_transplant, P = D_peers - C_elites, R = C_elites - C_opposite.
    Group differences use 95% percentile intervals from bootstrap resamples within each group.
    Returns (contrasts, per-cell means); both work on partial data and skip what is not yet computable."""
    import numpy as np
    import pandas as pd
    meta = {c["id"]: c for c in cases if "cell" in c}
    names = {"E": "elite_C_minus_B", "P": "peer_exposure", "R": "elite_direction_C_minus_opposite"}
    rec = per_case[per_case.metric == "recommendation"]
    steps = {}
    for arm in arms:
        pivot = rec[rec.arm == arm].pivot(index="case_id", columns="contrast", values="delta")
        steps[arm] = pd.DataFrame({k: pivot.get(v, pd.Series(dtype=float)) for k, v in names.items()})
        steps[arm]["EP"] = steps[arm]["E"] - steps[arm]["P"]
    if {"talkie", "qwen1930"} <= set(steps):
        steps["talkie_minus_qwen"] = steps["talkie"] - steps["qwen1930"]

    def kind(i): return meta[i]["cell"].split("_")[0]
    def source(i): return meta[i]["cell"].split("_")[1]
    tests = [  # (name, step, grouping, high, low, domain filter)
        ("X_primary", "EP", kind, "moral", "practical", lambda i: True),
        ("X_without_telephone2", "EP", kind, "moral", "practical", lambda i: i != "telephone2"),
        ("E_moral_minus_practical", "E", kind, "moral", "practical", lambda i: True),
        ("P_moral_minus_practical", "P", kind, "moral", "practical", lambda i: True),
        ("R_moral_minus_practical", "R", kind, "moral", "practical", lambda i: True),
        ("E_real_minus_fictional", "E", source, "real", "fictional", lambda i: True),
        ("E_known_failed_minus_not", "E", lambda i: meta[i].get("known_failed"), 1, 0,
         lambda i: source(i) == "fictional"),
        ("E_pre1931_minus_not", "E", lambda i: meta[i].get("pre1931"), 1, 0,
         lambda i: source(i) == "fictional"),
    ]
    rows = []
    for arm, frame in steps.items():
        for name, step, group, high, low, keep in tests:
            if arm == "talkie_minus_qwen" and name != "X_primary":
                continue
            z = frame[step].dropna()
            z = z[[i in meta and keep(i) for i in z.index]]
            hi = z[[group(i) == high for i in z.index]].to_numpy(dtype=float)
            lo = z[[group(i) == low for i in z.index]].to_numpy(dtype=float)
            if not len(hi) or not len(lo):
                continue
            draws = (rng.choice(hi, size=(n_boot, len(hi))).mean(axis=1)
                     - rng.choice(lo, size=(n_boot, len(lo))).mean(axis=1))
            rows.append({"arm": arm, "test": name, "estimate": hi.mean() - lo.mean(),
                         "ci_low": np.percentile(draws, 2.5), "ci_high": np.percentile(draws, 97.5),
                         "mean_high": hi.mean(), "mean_low": lo.mean(),
                         "n_high": len(hi), "n_low": len(lo)})
    cell_rows = []
    for arm, frame in steps.items():
        for step in ("E", "P", "R", "EP"):
            z = frame[step].dropna()
            z = z[[i in meta for i in z.index]]
            for cell, g in z.groupby(lambda i: meta[i]["cell"]):
                v = g.to_numpy(dtype=float)
                draws = rng.choice(v, size=(n_boot, len(v))).mean(axis=1)
                cell_rows.append({"arm": arm, "step": step, "cell": cell, "mean": v.mean(),
                                  "ci_low": np.percentile(draws, 2.5) if len(v) > 1 else np.nan,
                                  "ci_high": np.percentile(draws, 97.5) if len(v) > 1 else np.nan,
                                  "n_domains": len(v)})
    return pd.DataFrame(rows), pd.DataFrame(cell_rows)


def selftest():
    import tempfile
    from types import SimpleNamespace
    import pandas as pd
    import numpy as np
    a = ordered_trials(["talkie", "qwen1930"], EXTENDED, 2, 1930, CASES_BATTERY_A)
    assert len(a) == 8 * 9 * 2 * 2 == len({x["trial_id"] for x in a})
    assert a == ordered_trials(["talkie", "qwen1930"], EXTENDED, 2, 1930, CASES_BATTERY_A)
    for c in CASES:
        for condition in EXTENDED:
            v = vignette(c, condition)
            assert c["old"] in v if condition == "A_1930" else c["new"] in v
        assert c["explanation"] not in vignette(c, "C_elites")
        assert c["explanation"] in vignette(c, "C_explained")
        assert c["institution"] in vignette(c, "D_full")
        assert c["institution"] not in vignette(c, "D_media")
    sample = [{"turn":i+1,"question": QUESTIONS[i], "text": t} for i,t in enumerate([
        "I would choose the second course, as it is practical.",
        "I do not think either person is morally inferior.",
        "I would still choose the first course."])]
    judgment = {"recommendation":"B", "spontaneous_moral":"not_stated",
                "elicited_moral":"equal", "elicited_strength":0,
                "equalized_preference":"A", "reason_codes":["practical"],
                "evidence":{"recommendation":"I would choose the second course",
                            "spontaneous_moral":"",
                            "elicited_moral":"I do not think either person is morally inferior",
                            "equalized_preference":"I would still choose the first course"}}
    assert normalized_judgment(json.dumps(judgment), sample)["recommendation"] == "B"
    assert make_judge_payload({"case_id":"marriage","reverse":True,"answers":sample})["presented_first"] == "B"
    assert orientation("not_stated") != orientation("not_stated")
    # A complete synthetic fixture validates the analysis pipeline without claiming empirical findings.
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        cfg = dict(version=VERSION, batteries="a", cases=CASES_BATTERY_A,
                   arms=["talkie","qwen1930"], conditions=list(EXTENDED),
                   replicates=2, seed=1930, turns=3)
        (out / "manifest.json").write_text(json.dumps(cfg))
        effects = {"A_1930":0,"B_transplant":1,"C_elites":2,"D_full":5,
                   "C_explained":1,"C_opposite":0,"C_moral":4,"D_peers":3,"D_media":4}
        for tr in a:
            bump = effects[tr["condition"]]
            rec = "B" if bump >= 2 else "A"
            append_jsonl(out / "interviews.jsonl", {**tr,"status":"ok","answers":sample})
            append_jsonl(out / "judgments.jsonl", {"trial_id":tr["trial_id"],"status":"ok",
                "judgment": {"recommendation":rec,"spontaneous_moral":rec,
                             "elicited_moral":rec,"equalized_preference":rec,
                             "elicited_strength": int(bump >= 2),"reason_codes":["practical"],
                             "evidence":{}}})
        analyze(SimpleNamespace(out=tmp, bootstrap=100))
        for name in ("coded_observations.csv", "effects.csv", "case_effects.csv",
                     "model_delta_comparison.csv", "residual_ratios.csv",
                     "effects_recommendation.png"):
            assert (out / name).exists(), name
        df = pd.read_csv(out / "effects.csv")
        assert set(["ecology_B_minus_A", "elite_C_minus_B", "social_D_minus_C"]) <= set(df.contrast)
        assert np.isfinite(df.mean_delta).all()
    # Battery E: rated wording, selection and codes match battery_e/; five-world treatments; contrasts.
    import csv
    here = Path(__file__).resolve().parent / "battery_e"
    if here.exists():
        rated = {r["id"]: r for f in ("candidates.csv", "candidates_round2.csv")
                 for r in csv.DictReader(open(here / f, encoding="utf-8"))}
        chosen = {r["id"]: r["cell"] for r in csv.DictReader(open(here / "selection_round2.csv"))
                  if r["selected"] == "True"}
        codes = {r["id"]: r for r in csv.DictReader(open(here / "sources_round2.csv", encoding="utf-8"))}
        assert {c["id"]: c["cell"] for c in CASES_BATTERY_E} == chosen
        for c in CASES_BATTERY_E:
            assert all(c[k] == rated[c["id"]][k] for k in ("topic", "a", "b", "cell")), c["id"]
            if c["facts"] == "invented":
                assert (str(c["known_failed"]), str(c["pre1931"])) == (
                    codes[c["id"]]["known_failed_by_1930"], codes[c["id"]]["pre1931"]), c["id"]
        ext = {r["id"]: r for r in csv.DictReader(open(here / "extension_sentences.csv", encoding="utf-8"))}
        for c in CASES_BATTERY_E:  # deviation 4: approved wording for the four extra worlds
            assert (c["explanation"], c["institution"], c["cell"]) == (
                ext[c["id"]]["explanation"], ext[c["id"]]["institution"], ext[c["id"]]["cell"]), c["id"]
            assert ext[c["id"]]["status"] == "approved"
    ids = [c["id"] for c in CASES + CASES_BATTERY_E]
    assert len(ids) == len(set(ids))
    assert len({c["cell"] for c in CASES_BATTERY_E}) == 4 and len(CASES_BATTERY_E) == 40
    for c in CASES_BATTERY_E:
        assert c["facts"] == ("sourced" if c["cell"].endswith("real") else "invented")
        assert c["old"] in vignette(c, "A_1930") and c["new"] in vignette(c, "B_transplant")
        assert "There is one further fact" not in vignette(c, "B_transplant")
        assert c["a"] in exposure_text(c, "C_opposite") and c["b"] not in exposure_text(c, "C_opposite")
        assert exposure_text(c, "D_peers").startswith(exposure_text(c, "C_elites"))
        assert "practical reason is that " + c["explanation"] + ", and" in exposure_text(c, "C_explained")
        assert exposure_text(c, "D_full").startswith(exposure_text(c, "D_media"))
        assert exposure_text(c, "D_full").endswith(": " + c["institution"] + ".")
    rest = [w for w in EXTENDED if w not in FIVE]  # deviation 4: same presentation order as the five worlds
    order = {(x["case_id"], x["rep"]): x["reverse"] for x in ordered_trials(["talkie"], FIVE, 3, 1930, CASES_BATTERY_E)}
    r4 = ordered_trials(["talkie", "qwen1930"], rest, 3, 1930, CASES_BATTERY_E, "blocked")
    assert len(r4) == 960 and all(order[(x["case_id"], x["rep"])] == x["reverse"] for x in r4)
    e = ordered_trials(["talkie", "qwen1930"], FIVE, 3, 1930, CASES_BATTERY_E)
    assert len(e) == 1200 == len({x["trial_id"] for x in e})
    eb = ordered_trials(["talkie", "qwen1930"], FIVE, 3, 1930, CASES_BATTERY_E, "blocked")
    assert sorted(x["trial_id"] for x in eb) == sorted(x["trial_id"] for x in e)
    first = [x for x in eb if x["arm"] == "talkie"][:20]  # one domain per cell, all five worlds, rep 0
    assert {x["rep"] for x in first} == {0} and len({x["case_id"] for x in first}) == 4
    assert {get_case(x["case_id"])["cell"] for x in first} == {c["cell"] for c in CASES_BATTERY_E}
    assert {x["rep"] for x in [x for x in eb if x["arm"] == "talkie"][:200]} == {0}
    fake = []  # elite effect for moral domains, peer effect for practical ones: X = 1 - (-1) = 2
    for c in CASES_BATTERY_E:
        moral = c["cell"].startswith("moral")
        for arm in ("talkie", "qwen1930"):
            for contrast, d in (("elite_C_minus_B", 1 if moral else 0), ("peer_exposure", 0 if moral else 1),
                                ("elite_direction_C_minus_opposite", 1)):
                fake.append({"arm": arm, "contrast": contrast, "metric": "recommendation",
                             "case_id": c["id"], "delta": d, "n_pairs": 3})
    res, cells = battery_e_contrasts(pd.DataFrame(fake), CASES_BATTERY_E, ["talkie", "qwen1930"],
                                     np.random.default_rng(0), n_boot=200)
    res = res.set_index(["arm", "test"])
    assert cells.set_index(["arm", "step", "cell"]).loc[("talkie", "EP", "practical_real"), "mean"] == -1
    part, _ = battery_e_contrasts(pd.DataFrame(fake[:6]), CASES_BATTERY_E, ["talkie"], np.random.default_rng(0), 50)
    assert part.empty  # one domain: nothing to contrast yet, but no error
    assert res.loc[("talkie", "X_primary"), "estimate"] == 2
    assert res.loc[("talkie", "X_without_telephone2"), "n_low"] == 19
    assert res.loc[("talkie_minus_qwen", "X_primary"), "estimate"] == 0
    assert res.loc[("talkie", "R_moral_minus_practical"), "estimate"] == 0
    print("SELFTEST PASS: trial pairing, treatment contents, judge validation, synthetic analysis, battery E")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=["plan", "run", "judge", "analyze", "merge", "selftest"])
    p.add_argument("--out", default="hanson_results")
    p.add_argument("--sources", nargs="+", default=[], help="merge: run directories to combine into --out")
    p.add_argument("--design", choices=["core", "extended", "five"], default="extended")
    p.add_argument("--batteries", default="abc",
                   help="Case batteries to include: a (eco→B), b (more eco→B), c (eco→A drift test), d (Hanson-list changes, real facts), e (moral/practical x real/fictional; use --design five)")
    p.add_argument("--order", choices=["shuffled", "blocked"], default="shuffled",
                   help="blocked: complete (domain, replicate) blocks first, for early partial analysis")
    p.add_argument("--conditions", default="", help="Override with comma-separated condition IDs")
    p.add_argument("--arms", default="talkie,qwen1930")
    p.add_argument("--replicates", type=int, default=3)
    p.add_argument("--turns", type=int, default=3)
    p.add_argument("--seed", type=int, default=1930)
    p.add_argument("--local-url", default="http://127.0.0.1:8080/v1")
    p.add_argument("--local-model", default="")
    p.add_argument("--remote-model", default=DEFAULT_QWEN)
    p.add_argument("--judge-model", default=DEFAULT_JUDGE)
    p.add_argument("--temperature", type=float, default=0.5)
    p.add_argument("--max-tokens", type=int, default=430)
    p.add_argument("--local-max-tokens", type=int, default=180,
                   help="Per-turn cap for local Talkie (fits 2048 ctx with 3 turns)")
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--pause", type=float, default=0.15)
    p.add_argument("--limit", type=int, default=0, help="For pilot: max pending interviews/judgments")
    args = p.parse_args()
    if args.command == "selftest":
        selftest()
    elif args.command == "plan":
        cfg = get_config(args)
        total = len(cfg["cases"])*len(cfg["conditions"])*cfg["replicates"]*len(cfg["arms"])
        print(f"{total} interviews / up to {total*cfg['turns']} respondent completions + {total} judge calls")
        print("batteries:", cfg["batteries"], "domains:", len(cfg["cases"]))
        print("conditions:", ", ".join(cfg["conditions"]))
        print("sample vignette and 3 questions:\n")
        print(vignette(cfg["cases"][0], "B_transplant"))
        print("\n".join("Question: "+q for q in QUESTIONS))
    elif args.command == "run":
        run(args)
    elif args.command == "judge":
        judge(args)
    elif args.command == "analyze":
        analyze(args)
    elif args.command == "merge":
        merge(args)


if __name__ == "__main__":
    main()
