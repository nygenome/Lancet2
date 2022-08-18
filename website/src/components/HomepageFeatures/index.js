import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

const FeatureList = [
  {
    title: 'Localized',
    description: (
      <>
        With Lancet&apos;s localized assembly and construction of deBruijn graphs,
        Lancet is able to quickly and accurately detect variants in a tumor-normal 
        pair while working efficiently using as many CPU resources as available.
      </>
    ),
  },
  {
    title: 'User Friendly',
    description: (
      <>
        Standard Lancet variant calling only requires a tumor and normal sample along with
        an accompanying reference fasta and a designated path to output the vcf file to. 
        Check out the command line section for different options to customize a run.
      </>
    ),
  },
  {
    title: 'Joint Assembly',
    description: (
      <>
        As a tool, Lancet is able to jointly assemble reads from a tumor sample along with a 
        matched normal sample. This results in an increase in accuracy of identifying mutations, 
        especially indels, private to the tumor.
      </>
    ),
  },
];

function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center padding-horiz--md">
        <h3>{title}</h3>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
        <div className="col">
          <div className="text--center">
            <h2>Funding</h2>
            <p>Informatics Technology for Cancer Research (<a href="https://itcr.cancer.gov">ITCR</a>) under the NCI U01 award <a href="https://reporter.nih.gov/project-details/10304730">1U01CA253405-01A1</a>.</p>
          </div>
        </div>
      </div>
    </section>
  );
}
