// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

const lightCodeTheme = require('prism-react-renderer/themes/github');
const darkCodeTheme = require('prism-react-renderer/themes/dracula');

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'Lancet',
  staticDirectories: ['static'],
  tagline: 'Somatic variant caller with localized micro-assembly',
  url: 'https://nygenome.github.io',
  baseUrl: '/Lancet2/',
  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'img/favicon.ico',

  // GitHub pages deployment config.
  // If you aren't using GitHub pages, you don't need these.
  organizationName: 'nygenome', // Usually your GitHub org/user name.
  projectName: 'Lancet2', // Usually your repo name.

  // Even if you don't use internalization, you can use this field to set useful
  // metadata like html lang. For example, if your site is Chinese, you may want
  // to replace "en" with "zh-Hans".
  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          // Remove this to remove the "edit this page" links.
          editUrl:
            'https://github.com/nygenome/Lancet2/tree/main/website',
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          // Remove this to remove the "edit this page" links.
          editUrl:
            'https://github.com/nygenome/Lancet2/tree/main/website',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      }),
    ],
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      navbar: {
        title: 'Lancet',
        logo: {
          alt: 'Lancet Site Logo',
          src: 'img/nygc_logo.png',
        },
        items: [
          {
            type: 'doc',
            docId: 'getting_started',
            position: 'left',
            label: 'Docs',
          },
          {
            type: 'doc',
            docId: 'cli',
            position: 'left',
            label: 'Command Line',
          },
          {
            type: 'doc',
            docId: 'publications',
            position: 'left',
            label: 'Publications',
          },
          {to: '/blog', label: 'Blog', position: 'left'},
          {
            href: 'https://github.com/nygenome/Lancet2',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
      footer: {
        style: 'dark',
        links: [
          {
            title: 'Docs',
            items: [
              {
                label: 'Getting Started',
                to: '/docs/getting_started',
              },
              {
                label: 'Command line',
                to: '/docs/cli',
              },
            ],
          },
          {
            title: 'Community',
            items: [
              {
                label: 'GitHub',
                href: 'https://github.com/nygenome/Lancet2',
              },
              {
                label: 'Biostars',
                href: 'https://www.biostars.org/post/search/?query=lancet',
              },
            ],
          },
          {
            title: 'More',
            items: [
              {
                label: 'Blog',
                to: '/blog',
              },
              {
                label: 'Acknowledgments',
                to: '/docs/acknowledgment',
              }
            ],
          },
        ],
        copyright: `Copyright © ${new Date().getFullYear()} Lancet. Built with Docusaurus.`,
      },
      prism: {
        theme: lightCodeTheme,
        darkTheme: darkCodeTheme,
      },
    }),
};

module.exports = config;
