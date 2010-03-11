package edu.washington.gs.noble.crux.gui;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.logging.Logger;

import javax.swing.JPanel;

@SuppressWarnings("serial")
class ImagePanel extends JPanel {
	private BufferedImage image;
	private static Logger logger = Logger
	.getLogger("edu.washington.gs.noble.crux.gui");

	public ImagePanel(
			CruxComponentButton createIndex,
			CruxComponentButton search,
			CruxComponentButton computeQvalues,
			CruxComponentButton percolator,
			CruxComponentButton qranker
	) {
		super();
		// read the image
		java.net.URL imgURL = getClass().getResource("/edu/washington/gs/noble/crux/gui/schematic.png");
		try {
		  image = javax.imageio.ImageIO.read(imgURL);
		}
		catch (IOException exception) {
		  logger.info(exception.getMessage());
		  System.exit(-1);
		}
		
		setSize(new Dimension(image.getWidth(), image.getHeight()));
		setMaximumSize(new Dimension(image.getWidth(), image.getHeight()));
		setMinimumSize(new Dimension(image.getWidth(), image.getHeight()));
		setPreferredSize(new Dimension(image.getWidth(), image.getHeight()));
		setLayout(null);
		createIndex.setBounds(new Rectangle(88, 61, 170, 37));
		add(createIndex);
		search.setBounds(new Rectangle(400, 61, 170, 37));
		add(search);
		computeQvalues.setBounds(new Rectangle(790, 0, 160, 37));
		add(computeQvalues);
		percolator.setBounds(new Rectangle(816, 61, 120, 37));
		add(percolator);
		qranker.setBounds(new Rectangle(808, 126, 120, 37));
		add(qranker);
		/*ButtonGroup cruxButtons = new ButtonGroup();
		cruxButtons.add(createIndex);
		cruxButtons.add(search);
		cruxButtons.add(computeQvalues);
		cruxButtons.add(percolator);
		cruxButtons.add(qranker);
		cruxButtons.add(dummy);
		dummy.setSelected(true);*/
	}

	public void paintComponent(Graphics g) {
		if (image == null) {
			return;
		}
		g.drawImage(image, 0, 0, getWidth(), getHeight(), this);
	}
}